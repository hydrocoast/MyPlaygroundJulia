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
Δt = 1.0
t = 0.0:Δt:300.0
nt = length(t)

## phase
ϕ₁ = 0.0
ϕ₂ = 1/3*π

## wave parameters
a = 0.5
h = 10.0
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


## plot
anim = @animate for it = vcat(ones(Int64,10,1),collect(Int64,1:nt),nt*ones(Int64,10,1))
    @printf("%04d / %04d\r", it, nt)
    plt = plot(; label="", xlabel=L"x", ylabel=L"\eta", guidefont=12, tickfont=("Times",10), 
               title=@sprintf("%03d s",t[it]), titlefont=16, xlims=extrema(x), ylims=(-2a,2*a),
               size=(600,300), margin=2*Plots.PlotMeasures.mm)
    plt = plot!(plt, x, η₁[:,it]; label=L"\eta_1(k_1)", legendfontsize=10)
    plt = plot!(plt, x, η₂[:,it]; label=L"\eta_2(k_2=k_1+\Delta k)")
    plt = plot!(plt, x, η₁₂[:,it]; label=L"\eta_{1}+\eta_{2}")
    plt = plot!(plt, x, ηg[:,it]; lc=:gray, lw=1.0, label="")
    plt = plot!(plt, x, -ηg[:,it]; lc=:gray, lw=1.0, label="")
end

## save
gifname = "test_groupwaves.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=25)
