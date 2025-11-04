# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRIBase#fix-sim-block-issue
# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRICore#fix-sim-block-issue
using GLMakie
using KomaMRICore
using JLD2

## Params to get a slice profile of Δz = 6mm
zmax = 60e-3
Nspins = 20 # Use an even number to not have z = 0 in the grid
z = range(-zmax, zmax, Nspins) |> collect
xs = copy(z)
ys = copy(z)
x = [x for (x, y) in Iterators.product(xs, ys)][:]
y = [y for (x, y) in Iterators.product(xs, ys)][:]
x_grid_mm = reshape(x .* 1e3, length(xs), length(ys))
y_grid_mm = reshape(y .* 1e3, length(xs), length(ys))
heatmap_plane_z = fill(-zmax * 1e3, size(x_grid_mm))
plane_shape = size(x_grid_mm)

# Setting up input objects
sys = Scanner()
seq = load(joinpath(@__DIR__, "JuliaLogo_waveform_resampled.jld2"))["seq_resampled"]
seq.ADC[1] = ADC(120, dur(seq)/2, dur(seq)/2)
seq = (3.0 .+ 0im) * seq # Scale B1 amplitude
obj = Phantom(; x, y)
# plot_seq(seq; show_adc=true)

## Simulate
sim_params_gt = KomaMRICore.default_sim_params()
sim_params_gt["return_type"] = "mat"
sim_params_gt["sim_method"] = BlochDict(save_Mz=true)
sim_params_gt["Nthreads"] = 1 # Problems with BlochDict and multithreading for !=1
mag = simulate(obj, seq, sys; sim_params=sim_params_gt)

## Plot
i = Observable(1)
us = @lift(real.(mag[$i, :, 1]))
vs = @lift(imag.(mag[$i, :, 1]))
ws = @lift(real.(mag[$i, :, 2]))
saturation = @lift(abs.(mag[$i, :, 1]))
saturation2 = @lift(-3 .* abs.(mag[$i, :, 1]))
hue = @lift(angle.(mag[$i, :, 1]))
seqd = discretize(seq)
t = seqd.t .* 1e3
t_adc = t[seqd.ADC]
current_time = @lift(t_adc[$i])
B1 = seqd.B1 .* 1e6 .* 0.8
Gx = seqd.Gx .* 1e3
Gy = seqd.Gy .* 1e3
Gz = seqd.Gz .* 1e3

# Show magnetic field 
# Gx_adc = Gx[seqd.ADC]
# Gy_adc = Gy[seqd.ADC]
# Gz_adc = Gz[seqd.ADC]
# B_plane = @lift(reshape(abs.(x .* Gx_adc[$i] .+ y .* Gy_adc[$i]).^2, plane_shape...))

fig = Figure(size = (1800, 1800), fontsize=40, backgroundcolor = :black)
ax = Axis3(fig[2, 1], aspect = (1, 1, 1), xlabel=L"x", ylabel=L"y", zlabel=L"z", backgroundcolor=:black)
xlims!(ax, -zmax * 1e3, zmax * 1e3)
ylims!(ax, -zmax * 1e3, zmax * 1e3)
zlims!(ax, -zmax * 1e3, zmax * 1e3)
# ax.yreversed = true
for k in 1:4
    setproperty!(ax, Symbol("xspinecolor_$k"), :white)
    setproperty!(ax, Symbol("yspinecolor_$k"), :white)
    setproperty!(ax, Symbol("zspinecolor_$k"), :white)
end
ax.xlabelcolor = :white
ax.ylabelcolor = :white
ax.zlabelcolor = :white
ax.xlabelsize = 90
ax.ylabelsize = 90
ax.zlabelsize = 90

# Show magnetic field plane
# surface!(ax, x_grid_mm, y_grid_mm, heatmap_plane_z; color=B_plane, colormap=:grays, shading=false)

ar = arrows3d!(ax,
    x .* 1e3,
    y .* 1e3,
    zeros(length(x)) .* 1e3,
    us, vs, ws;
    lengthscale = 6,
    color = saturation,
    align = :center,
    colorrange=(0.0, 1.0),
    # colormap=:grays
)
# lines!(ax, 1e3 .* zmax .* ones(size(z)), saturation2, z .* 1e3; color=saturation, colorrange=(0.0, 1.0), linewidth=5)

ax2 = Axis(fig[1, 1], backgroundcolor = :black)
vlines!(ax2, current_time, color=:white, linewidth=5)
lines!(ax2, t[B1 .!= 0], real.(B1[B1 .!= 0]), linewidth=5, label=L"B_1")
lines!(ax2, t, Gx, linewidth=5, label=L"G_x")
lines!(ax2, t, Gy, linewidth=5, label=L"G_y")
xlims!(ax2, maximum(t)/2, maximum(t))
hidespines!(ax2)
hidedecorations!(ax2)
leg = axislegend(ax2, labelcolor=:white, backgroundcolor=(:black, 1), framecolor=:black, labelsize=90)
fig[1, 1] = ax2
fig[1, 2] = leg
# hidedecorations!(ax)
rowsize!(fig.layout, 1, Auto(1))
rowsize!(fig.layout, 2, Auto(5))
display(fig)

## Animate
nframes = size(mag, 1)
frames0 = collect(1:nframes)
frames = [frames0; last(frames0) * ones(40)] # Hold the last frame for a while
gif_path = joinpath(@__DIR__, "figures", "rf_2dexcitation_logo.gif")
record(fig, gif_path, frames; framerate = 20) do k
    i[] = k
end

## Interactive
sl_x = Slider(fig[3, 1], range = frames0, startvalue = 1, linewidth=40)
lift(sl_x.value) do k
    i[] = Int(k)
end
fig
