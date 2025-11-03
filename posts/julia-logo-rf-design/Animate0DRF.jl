# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRIBase#fix-sim-block-issue
# add https://github.com/JuliaHealth/KomaMRI.jl:KomaMRICore#fix-sim-block-issue
using GLMakie
using KomaMRICore

## Params to get a slice profile of Δz = 6mm
Trf = 3.2e-3
B1 = (π/2) / (2π*γ*Trf) * 2
Δf = 5 / Trf
# Setting up input objects
sys = Scanner()
seq = Sequence()
seq += RF(-B1 .* cos.(-2π * Δf * range(0, Trf, 400)), Trf)
seq.ADC[1] = ADC(120, dur(seq))
obj = Phantom(; x=[0.0], Δw=[2π * Δf])

## Simulate
sim_params_gt = KomaMRICore.default_sim_params()
sim_params_gt["return_type"] = "mat"
sim_params_gt["sim_method"] = BlochDict(save_Mz=true)
sim_params_gt["Nthreads"] = 1 # Problems with BlochDict and multithreading for !=1
mag = simulate(obj, seq, sys; sim_params=sim_params_gt)

## Plot
i = Observable(1)
seqd = discretize(seq)
us = @lift(real.(mag[$i, :, 1]))
vs = @lift(imag.(mag[$i, :, 1]))
ws = @lift(real.(mag[$i, :, 2]))
saturation = @lift(abs.(mag[$i, :, 1]))
saturation2 = @lift(3 .* abs.(mag[$i, :, 1]))
hue = @lift(angle.(mag[$i, :, 1]))

fig = Figure(size = (1800, 1800), fontsize=40, backgroundcolor = :black)
ax = Axis3(fig[2, 1], aspect = (1, 1, 1), xlabel=L"x", ylabel=L"y", zlabel=L"z", backgroundcolor=:black)
ar = arrows3d!(ax,
    obj.x .* 1e3,
    obj.y .* 1e3,
    obj.z .* 1e3,
    us, 
    vs, 
    ws;
    lengthscale = 0.4,
    color = saturation,
    align = :tail,
    colorrange=(0.0, 1.0),
)
B1_adc = real.(KomaMRIBase.get_rfs(seq, seqd.t[seqd.ADC])[1] * 0.2e6)
b1 = @lift([B1_adc[$i]])
ar = arrows3d!(ax,
    obj.x .* 1e3,
    obj.y .* 1e3,
    obj.z .* 1e3,
    b1,
    [0.0],
    [0.0];
    lengthscale = 0.4,
    color = Makie.wong_colors()[1],
    align = :tail,
    colorrange=(0.0, 1.0),
)
xlims!(ax, -0.5, 0.5)
ylims!(ax, -0.5, 0.5)
zlims!(ax, -0.5, 0.5)
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
# lines!(ax, 1e3 .* zmax .* ones(size(z)), saturation2, z .* 1e3; color=saturation, colorrange=(0.0, 1.0), linewidth=5)

ax2 = Axis(fig[1, 1], backgroundcolor = :black)
seqd = discretize(seq)
t = seqd.t .* 1e3
t_adc = t[seqd.ADC]
current_time = @lift(t_adc[$i])
B1 = seqd.B1 .* 1e6 .* 0.8
Gz = seqd.Gz .* 1e3
vlines!(ax2, current_time, color=:white, linewidth=5)
lines!(ax2, t, real.(B1), linewidth=5, label=L"B_1", color=Makie.wong_colors()[1])
# lines!(ax2, t, real.(B1), linewidth=5, label=L"B_{1,x}", linestyle=:dash, color=Makie.wong_colors()[1])
# lines!(ax2, t, Gz, linewidth=5, label=L"G_z")
ylims!(ax2, -1.2 * maximum(abs.(B1)), 1.2 * maximum(abs.(B1)))
hidespines!(ax2)
hidedecorations!(ax2)
leg = axislegend(ax2, labelcolor=:white, backgroundcolor=(:black, 1), framecolor=:black, labelsize=90)
fig[1, 1] = ax2
fig[1, 2] = leg
rowsize!(fig.layout, 1, Auto(1))
rowsize!(fig.layout, 2, Auto(5))
display(fig)

## Animate
nframes = size(mag, 1)
frames0 = collect(1:nframes)
frames = [frames0; last(frames0) * ones(40)] # Hold the last frame for a while
record(fig, "figures/rf_excitation_onespin.gif", frames; framerate = 20) do k
    i[] = k
end

## Interactive
nframes = size(mag, 1)
frames0 = collect(1:nframes)
sl_x = Slider(fig[3, 1], range = frames0, startvalue = 1, linewidth=40)
lift(sl_x.value) do k
    i[] = Int(k)
end
fig