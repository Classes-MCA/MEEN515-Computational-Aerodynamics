# include("Desktop/School/Aerodynamics/Homework 2/Homework 2-1.jl")

using Pkg
Pkg.add("PyPlot")
using PyPlot

include("taf.jl")

e = 2/100
p = 4/10
t = 12/100

c = 1
alpha = 2*pi/180 # Example Fourier Coefficients based on AOA = 2 Degrees
dphi = 0.001

phi = 0:dphi:pi
x = zeros(length(phi),1)
for i = 1:length(phi)
    x[i] = c/2 * (1-cos(phi[i]))
end

T, ybar, dTdx, dybardx = naca4(e, p, t, x)

# Plot the airfoil
figure()
plot(x,ybar+0.5*T)
plot(x,ybar-0.5*T)
ylim(-.3,.3)

integrand = zeros(1,length(phi)) # to hold the integral values before they are summed

# Solving for A0
for i = 1:length(phi)
    integrand[i] = dybardx[i]*dphi
end

A0 = alpha - 1/pi * sum(integrand)

############################################
# Solving for the other Fourier Coefficients
coefficients = zeros(1,4)

for n = 1:length(coefficients)
    for i = 1:length(phi)
        integrand[i] = dybardx[i] * cos(n * phi[i]) * dphi
    end

    coefficients[n] = 2/pi * sum(integrand)

end

println("A0 = ",A0) # .0304
println("A1 = ",coefficients[1]) # 0.08149
println("A2 = ",coefficients[2]) # 0.01386
println("A3 = ",coefficients[3]) # 0.00277
println("A4 = ",coefficients[4]) # -0.0021