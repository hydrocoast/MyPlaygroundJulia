using Plots
using Printf
using LaTeXStrings

## constant
g = 9.8

## spatial 
Δx = 1/64*π
x = -10π:Δx:10π
nx = length(x)

## temporal
Δt = 0.25
t = 0.0:Δt:100.0
nt = length(t)

## phase
ϕ₁ = 0.0
ϕ₂ = 1/3*π

## wave parameters
a = 0.5
h = 20.0
k₁ = 5.0
k₂ = k₁ + 0.1

## dispersion relation
ω(k,h) = sqrt(g*k*tanh(k*h))
ω₁ = ω(k₁,h)
ω₂ = ω(k₂,h)
T₁ = 2π/ω₁
T₂ = 2π/ω₂
f₁ = 1/T₁
f₂ = 1/T₂
c₁ = ω₁/k₁
c₂ = ω₂/k₂
cp = (ω₁+ω₂)/(k₁+k₂)
cg = abs(ω₁-ω₂)/abs(k₁-k₂)

## preallocate
η₁ = zeros(nx,nt)
η₂ = similar(η₁)
η₁₂ = similar(η₁)
ηg = similar(η₁)


## calc
for j = 1:nt
for i = 1:nx
    η₁[i,j] = a*sin(k₁*x[i] - ω₁*t[j] + ϕ₁)
    η₂[i,j] = a*sin(k₂*x[i] - ω₂*t[j] + ϕ₂)
    η₁₂[i,j] = η₁[i,j] + η₂[i,j]
    ηg[i,j] = 2*a*cos(0.5*((k₁-k₂)*x[i] - (ω₁-ω₂)*t[j] + (ϕ₁-ϕ₂)))
end
end

_, ind1 = findmax(η₁[:,1], dims=1)
_, ind2 = findmax(η₂[1:ind1[1]+20,1], dims=1)
_, ind12 = findmax(η₁₂[:,1], dims=1)


c = palette(:auto)
## plot
anim = @animate for it = vcat(ones(Int64,10,1),collect(Int64,1:nt),nt*ones(Int64,10,1))
#anim = @animate for it = vcat(ones(Int64,5,1),collect(Int64,1:100),100*ones(Int64,5,1))
    @printf("%04d / %04d\r", it, nt)
    plt = plot(; label="", xlabel=L"x", ylabel=L"\eta", guidefont=12, tickfont=("Times",10), 
               title=@sprintf("%05.1f s",t[it]), titlefont=14, xlims=extrema(x), ylims=(-2a-0.05,2*a+0.05),
               size=(800,400), margin=3*Plots.PlotMeasures.mm)
    plt = plot!(plt, x, η₁[:,it]; label=L"\eta_1(k_1)", legendfontsize=10, legend=:bottomright)
    plt = plot!(plt, x, η₂[:,it]; label=L"\eta_2(k_2=k_1+\Delta k)")
    plt = plot!(plt, x, η₁₂[:,it]; label=L"\eta_{1}+\eta_{2}")
    plt = plot!(plt, x, ηg[:,it]; lc=:gray, lw=1.0, label="")
    plt = plot!(plt, x, -ηg[:,it]; lc=:gray, lw=1.0, label="")

    # celerity marker
    xm1 = x[ind1[1]]+c₁*t[it]
    while xm1 > x[end]; xm1 = xm1 - (x[end]-x[1]); end
    xm2 = x[ind2[1]]+c₂*t[it]
    while xm2 > x[end]; xm2 = xm2 - (x[end]-x[1]); end
    xmg = x[ind12[1]]+cg*t[it]
    while xmg > x[end]; xmg = xmg - (x[end]-x[1]); end
    xmp = x[ind12[1]]+cp*t[it]
    while xmp > x[end]; xmp = xmp - (x[end]-x[1]); end
    _, ind_xmp = findmin(abs.(x.-xmp), dims=1)
    plt = plot!(plt, [xm1], [a], m=:circle, color=c[1], label="", markersize=6)
    plt = plot!(plt, [xm2], [a], m=:circle, color=c[2], label="", markersize=6)
    plt = plot!(plt, [xmg], [2a], m=:circle, color=:gray, label="", markersize=6)
    plt = plot!(plt, [x[ind_xmp[1]]], [η₁₂[ind_xmp,it]], m=:circle, color=c[3], label="", markersize=6)
end

## save
gifname = "test_groupvelocity.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=20)
