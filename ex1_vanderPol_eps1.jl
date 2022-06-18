using Plots
using LaTeXStrings
using Printf

ll = L" \ddot{x} - \varepsilon (1-x^2) \dot{x} + x = 0 ~~~ (\varepsilon > 0)"

ε = 0.1
x₀ = 0.01
v₀ = 0.0
a₀ = 0.0

t₀ = 0.0
Δt = 0.01
tend = 150.0
time_sim = t₀:Δt:tend
nt = length(time_sim)


x = zeros(nt,1)
v = zeros(nt,1)
a = zeros(nt,1)
x[1] = x₀
v[1] = v₀
a[1] = a₀
for it = 2:nt
    x[it] = x[it-1] + v[it-1]*Δt
    v[it] = v[it-1] + a[it-1]*Δt
    a[it] = ε*(1-x[it-1]^2)*v[it-1] - x[it-1]
end

#=
plt1 = plot(time_sim, x; label="", xlabel=L"t", ylabel=L"x", guidefont=16, tickfont=("Times",12), title=ll, titlefont=16)
plt1 = annotate!(plt1, 35, 1.8, Plots.text(@sprintf("\$ \\varepsilon = %0.2f\$",ε), :right, :middle))
plt1 = annotate!(plt1, 35, 1.4, Plots.text(@sprintf("\$ x _0 = %0.2f\$",x₀), :right, :middle))
plt2 = plot(x,v; label="", xlabel=L"x", ylabel=L"v", guidefont=16, tickfont=("Times",12), title="Van der Pol oscillator", titlefont=("Times",16))
plt = plot(plt2, plt1, layout=(1,2), size=(1200,600))
=#


anim = @animate for it = vcat(ones(Int64,4,1),collect(Int64,1:100:nt),nt*ones(Int64,10,1))
    @printf("%05d / %05d\r", it, nt)
    pltA = plot(time_sim[1:it], x[1:it]; label="", xlabel=L"t", ylabel=L"x", guidefont=16, tickfont=("Times",12), title=ll, titlefont=16)
    pltA = plot!(pltA, [time_sim[it]], [x[it]]; label="", marker=:circle, mfc=:orange, ms=6, xlims=(-2,152), ylims=(-2.2,2.2))
    pltA = annotate!(pltA, 35, 1.8, Plots.text(@sprintf("\$ \\varepsilon = %0.2f\$",ε), :right, :middle))
    pltA = annotate!(pltA, 35, 1.4, Plots.text(@sprintf("\$ x _0 = %0.2f\$",x₀), :right, :middle))
    pltB = plot(x[1:it],v[1:it]; label="", xlabel=L"x", ylabel=L"v", guidefont=16, tickfont=("Times",12), title="Van der Pol oscillator", titlefont=("Times",16))
    pltB = plot!(pltB, [x[it]],[v[it]]; label="", marker=:circle, mfc=:orange, ms=6, xlims=(-2.2,2.2), ylims=(-2.2,2.2))
    pltG = plot(pltB, pltA, layout=(1,2), size=(1200,600))
end

## save
gifname = "test_vanderPol_eps1.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=8)
