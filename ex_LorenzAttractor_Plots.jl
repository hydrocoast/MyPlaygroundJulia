#=
Original code
https://docs.juliaplots.org/latest/user_gallery/misc/gr_lorenz_attractor/#Lorenz-Attractor
=#

using Printf
using LaTeXStrings
using Plots
gr()

# define the Lorenz attractor
Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8 / 3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x)
    dy = l.x * (l.ρ - l.z) - l.y
    dz = l.x * l.y - l.β * l.z
    l.x += l.dt * dx
    l.y += l.dt * dy
    l.z += l.dt * dz
end

attractor = Lorenz()

#=
# initialize a 3D plot with 1 empty series
plt = plot3d(
    1,
    xlim = (-30, 30),
    ylim = (-30, 30),
    zlim = (0, 55),
    title = "Lorenz Attractor",
    marker = :none,
    label = "",
    view = (45,45),
    margin = 3*Plots.PlotMeasures.mm,
    size = (600,600),
    xlabel = L"x", ylabel=L"y", zlabel=L"z"
)
=#

# build an animated gif by pushing new points to the plot, saving every 10th frame
# equivalently, you can use `@gif` to replace `@animate` and thus no need to explicitly call `gif(anim)`.
gifname = "Lorenz.gif"

#=
anim = @animate for i = 1:1500
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
end every 10
=#

nt = 3001
x = zeros(Float64, nt, 1)
y = similar(x)
z = similar(x)
t = attractor.dt .* collect(0:nt-1)

for i = 1:nt
    step!(attractor)
    x[i], y[i], z[i] = attractor.x, attractor.y, attractor.z
end


anim = @animate for i = 1:nt
    @printf("%d/%d\r",i,nt)
    plt = plot3d(x[1:i], y[1:i], z[1:i], 
    xlim = (-30, 30), ylim = (-30, 30), zlim = (0, 55),
    title = @sprintf("Lorenz Attractor\n%4.1f",attractor.dt*(i-1)),
    titlefont = Plots.font("Monospace", 14),
    marker = :none,
    line_z = t, colorbar = false,
    size = (600,600),
    xlabel = L"x", ylabel=L"y", zlabel=L"z", labelfontsize = 16,
    tickfont = Plots.font("Monospace", 10),
    label = "", view = (20,90),
    margin = 3*Plots.PlotMeasures.mm,
    )

end every 20

gif(anim, gifname; fps=20)
