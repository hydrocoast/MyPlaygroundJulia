using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Printf

Ldxdt = L"\frac{dx}{dt} = \varepsilon_1 x - k_1xy"
Ldydt = L"\frac{dy}{dt} = -\varepsilon_2 x + k_2xy"

ε₁ = 8.0
k₁ = 4.0
ε₂ = 18.0
k₂ = 6.0

x₀ = 8
y₀ = 3

t₀ = 0.0
Δt = 0.001
tend = 3.0
time_sim = t₀:Δt:tend
nt = length(time_sim)


x = zeros(nt,1)
y = zeros(nt,1)
x[1] = x₀
y[1] = y₀
for it = 2:nt
    x[it] = x[it-1] + ( ε₁*x[it-1] - k₁*x[it-1]*y[it-1])*Δt
    y[it] = y[it-1] + (-ε₂*y[it-1] + k₂*x[it-1]*y[it-1])*Δt
end

c = palette(:auto)


plt1 = plot(time_sim, hcat(x,y); label="", xlabel=L"t", ylabel=L"x,y", guidefont=16, tickfont=("Times",12), titlefont=16, xlims=(-0.2,3.2), ylims=(-0.2,10.2))
plt1 = annotate!(plt1, 0.2, 10.0, Plots.text(Ldxdt, :left, :middle, c[1]))
plt1 = annotate!(plt1, 0.2, 9.00, Plots.text(Ldydt, :left, :middle, c[2]))
plt2 = plot(x,y; lc=:green, label="", xlabel=L"x", ylabel=L"y", guidefont=16, tickfont=("Times",12), title="Lotka-Volterra equations", titlefont=("Times",16), xlims=(-0.2,10.2), ylims=(-0.2,10.2))
plt = plot(plt2, plt1, layout=(1,2), size=(1200,600))


anim = @animate for it = vcat(ones(Int64,4,1),collect(Int64,1:50:nt),nt*ones(Int64,10,1))
    @printf("%05d / %05d\r", it, nt)
    plt1 = plot(time_sim[1:it], hcat(x[1:it],y[1:it]); label="", xlabel=L"t", ylabel=L"x,y", guidefont=16, tickfont=("Times",12), titlefont=16, xlims=(-0.2,3.2), ylims=(-0.2,10.2))
    plt1 = annotate!(plt1, 0.2, 10.0, Plots.text(Ldxdt, :left, :middle, c[1]))
    plt1 = annotate!(plt1, 0.2, 9.00, Plots.text(Ldydt, :left, :middle, c[2]))
    plt2 = plot(x[1:it],y[1:it]; lc=:green, label="", xlabel=L"x", ylabel=L"y", guidefont=16, tickfont=("Times",12), title="Lotka-Volterra equations", titlefont=("Times",16), xlims=(-0.2,10.2), ylims=(-0.2,10.2))
    plt2 = plot!(plt2, [x[it]],[y[it]]; label="", marker=:circle, mc=:green, ms=6)
    pltG = plot(plt2, plt1, layout=(1,2), size=(1200,600), left_margin=5Plots.PlotMeasures.mm, bottom_margin=5Plots.PlotMeasures.mm)
end

## save
gifname = "test_LotkaVolterra.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=8)
