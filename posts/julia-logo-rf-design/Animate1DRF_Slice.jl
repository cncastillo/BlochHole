# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRIBase#fix-sim-block-issue
# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRICore#fix-sim-block-issue
using GLMakie
using KomaMRICore

## Params to get a slice profile of Δz = 6mm
B1 = 4.9e-6 * 2 * 3 / 2
Trf = 3.2e-3
TBP = 8
Δz = 4e-3
zmax = 5e-3
fmax = TBP / Trf
Nspins = 40 # Use an even number to not have z = 0 in the grid
z = range(-zmax, zmax, Nspins) |> collect
Gz = fmax / (γ * Δz)
f = γ * Gz * z
xs = copy(z)
zs = copy(z)
x = [x for (x, z) in Iterators.product(xs, zs)][:]
z = [z for (x, z) in Iterators.product(xs, zs)][:]

# Setting up input objects
sys = Scanner()
seq = PulseDesigner.RF_sinc(B1 .+ 0im, Trf, sys; G=[0; 0; Gz], TBP=TBP)
seq.ADC[1] = ADC(120, dur(seq))
obj = Phantom(; x, z)
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

ar = arrows3d!(ax,
    x .* 1e3,
    zeros(length(x)) .* 1e3,
    z .* 1e3,
    us, vs, ws;
    lengthscale = 0.4,
    color = saturation,
    align = :center,
    colorrange=(0.0, 1.0),
    # colormap=:grays
)
lines!(ax, 1e3 .* zmax .* ones(size(z)), saturation2, z .* 1e3; color=saturation, colorrange=(0.0, 1.0), linewidth=5)

ax2 = Axis(fig[1, 1], backgroundcolor = :black)
seqd = discretize(seq)
t = seqd.t .* 1e3
t_adc = t[seqd.ADC]
current_time = @lift(t_adc[$i])
B1 = seqd.B1 .* 1e6 .* 0.8
Gz = seqd.Gz .* 1e3
vlines!(ax2, current_time, color=:white, linewidth=5)
lines!(ax2, t[B1 .!= 0], real.(B1[B1 .!= 0]), linewidth=5, label=L"B_1")
lines!(ax2, t, Gz, linewidth=5, label=L"G_z")
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
gif_path = joinpath(@__DIR__, "figures", "rf_excitation.gif")
record(fig, gif_path, frames; framerate = 20) do k
    i[] = k
end

## Interactive
sl_x = Slider(fig[3, 1], range = frames0, startvalue = 1, linewidth=40)
lift(sl_x.value) do k
    i[] = Int(k)
end
fig
