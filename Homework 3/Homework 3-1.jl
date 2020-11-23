# This program calculates the boundary layer details for a flat plate using two different
# methods and comparing them with standard solutions

# include("Box/Classes/Aerodynamics/Homework 3/Homework 3-1.jl")

using Pkg
Pkg.add("PyPlot")
#Pkg.add("DifferentialEquations")
using PyPlot
#using DifferentialEquations


Rex = 10^6; # Reynolds number based on length of the plate
nu = 20*10^-6; # This term shows up looking like a velocity in the equations, but it's a nu
xInitial = 10^-6; # Numerically we can't start at exactly x = 0, will result in Inf
dx = 10^(-5);
plateLength = 1; # meters
x = xInitial:dx:plateLength;
dVedx = zeros(length(x)); # For a flat plate, this is all zeros. Not true for a curved surface or a pressure gradient

edgeVelocity = Rex * nu / plateLength; # Will be needed in several places. For a flat plate this is the freestream velocity

#########################################################
# Method 1: Thwaite's Method - Laminar and Incompressible

# Numerically solve the ODE for theta (momentum thickness)
Rex0 = edgeVelocity * x[1] / nu
thetaInitial = 0.664 * x[1] / sqrt(Rex0); # FIXME: Should this Reynolds number be any different from Rex?

# Solving the ODE with a for loop
theta = zeros(length(x));
theta[1] = thetaInitial;
dThetadx = zeros(length(x))
for i = 1:length(x) - 1

    dThetadx[i] = 0.225 * nu / (edgeVelocity * theta[i]) - 3 * (theta[i] / edgeVelocity) * dVedx[i]

    theta[i + 1] = theta[i] + (dThetadx[i] * dx)

end

# Setting the last value of dThetadx
dThetadx[end] = 0.225 * nu / (edgeVelocity * theta[end]) - 3 * (theta[end] / edgeVelocity) * dVedx[end]

## Calculate lambda
lambda = theta.^2 ./ nu .* dVedx; # FIXME: All zeros? dVedx = [0...0] for a flat plate with no pressure gradient


## Calculate l and H
# These equations are empirical
l = zeros(length(x))
H = zeros(length(x))
for i = 1:length(lambda)
    if lambda[i] >= 0 && lambda[i] <= 0.1
        l[i] = 0.22 + 1.57*lambda[i] -1.8*lambda[i]^2
        H[i] = 2.61 - 3.75*lambda[i] -5.24*lambda[i]^2
    elseif lambda[i] > -0.1 && lambda[i] < 0
        l[i] = 0.22 + 1.402*lambda[i] + 0.018*lambda[i] / (0.107 + lambda[i])
        H[i] = 2.088 + 0.0731 / (0.14 + lambda[i])
    elseif lambda[i] > 0.1
        l[i] = 0.22 + 1.57*0.1 -1.8*0.1^2
        H[i] = 2.61 -3.75*0.1 -5.24*0.1^2
    elseif lambda[i] < -0.1
        println("Laminar Separation predicted at x = ", x[i])
    end
end

## Calculate delta star (displacement thickness)
displacementThickness = H.*theta

## Calculate cf (local skin friction coefficient)
cf = 2 .* (dThetadx .+ dVedx.*theta./edgeVelocity .* (H .+ 2))

## Calcualte Cf (total skin friction coefficient, integrate cf to find this)
Cf = sum(cf.*dx) # FIXME: This varied as a function of plate length while Blasius didn't
println("Thwaite's Method Cf = ",Cf)

## Storing values to plot for Thwaite's Method
theta_Thwaite = theta
displacementThickness_Thwaite = displacementThickness
H_Thwaite = H
cf_Thwaite = cf

#################################################
# Comparison: Blasius Solution

Rex = edgeVelocity .* x ./ nu
ReL = edgeVelocity * plateLength / nu
displacementThickness_Blasius = zeros(length(x))
theta_Blasius = zeros(length(x))
cf_Blasius = zeros(length(x))
for i = 1:length(x)
    displacementThickness_Blasius[i] = 1.72 * x[i] / sqrt(Rex[i])
    theta_Blasius[i] = 0.664 * x[i] / sqrt(Rex[i])
    cf_Blasius[i] = 0.664 / sqrt(Rex[i])
end

Cf = 1.328 / sqrt(ReL)
println("Blasius Solution Cf = ",Cf)

## Storing values to plot for Blasius Solution
theta_Blasius = theta_Blasius
displacementThickness_Blasius = displacementThickness_Blasius
H_Blasius = displacementThickness_Blasius./theta_Blasius
cf_Blasius = cf_Blasius

#################################################
# Head's Method for turbulent boundary layers

# We will be implementing the Julia DifferentialEquations.jl solver

# # Example differential equations
# function lorenz!(du,u,p,t)
#     du[1] = 10.0*(u[2]-u[1])
#     du[2] = u[1]*(28.0-u[3]) - u[2]
#     du[3] = u[1]*u[2] - (8/3)*u[3]
# end

# u0 = [1.0;0.0;0.0]
# tspan = (0.0,100.0)
# prob = ODEProblem(lorenz!,u0,tspan)
# sol = solve(prob)

# Define a function that represents the equations we'll be using
# function HeadsMethod!(du,u,p,t)

#     du[1] = 0.0306/u[2] * (u[1] - 3)^(-0.6169) - dVedx*u[1]/edgeVelocity - du[2]*u[1]/u[2]
    
#     #if u[1] > 5.
#     cf = 1; H = 4;
#     du[2] = cf/2 - dVedx*u[2]/edgeVelocity*(H + 2)

# end

# u0 = [1.0;0.0;0.0]
# xspan = (0.0,10.0)
# prob = ODEProblem(HeadsMethod!,u0,xspan)
# sol = solve(prob)

# Brute force Head's Method
# I will substitute in the equation for dthetadx into the equation for dH1dx

# Creating arrays and setting intial values
thetaInitial = 10^(-6)
HInitial     = 1.4
theta = zeros(length(x))
H1    = zeros(length(x))
H     = zeros(length(x))
cf    = zeros(length(x))
theta[1] = thetaInitial
H[1]     = HInitial

# Solving for the intial H1 value
if H[1] <= 1.6
    H1[1] = 0.8234*(H[1] - 1.1)^(-1.287) + 3.3
else
    H1[1] = 1.5501*(H[1] - 0.6778)^(-3.064) + 3.3
end

for i = 1:length(x)-1

    # Solving for H
    if H1[i] >= 5.3
        H[i] = 0.86*(H1[i] - 3.3)^(-0.777) + 1.1
    else
        H[i] = 1.1538*(H1[i] - 3.3)^(-0.326) + 0.6778
    end

    # Solving for cf
    cf[i] = 0.246*10^(-0.678*H[i]) * (edgeVelocity * theta[i] / nu)^(-0.268)

    # Solving for dThetadx
    dThetadx = cf[i]/2 - dVedx[i]*theta[i]/edgeVelocity * (H[i] + 2)

    # Solving for dH1dx
    dH1dx = 0.0306/theta[i] * (H1[i] - 3)^(-0.6169) - dVedx[i]*H1[i]/edgeVelocity - dThetadx*H1[i]/theta[i]

    theta[i + 1] = theta[i] + dThetadx*dx
    H1[i + 1]    = H1[i] + dH1dx*dx

end

# Solving for final H value
if H1[end] >= 5.3
    H[end] = 0.86*(H1[end] - 3.3)^(-0.777) + 1.1
else
    H[end] = 1.1538*(H1[end] - 3.3)^(-0.326) + 0.6778
end

# Solving for final cf value
cf[end] = 0.246*10^(-0.678*H[end]) * (edgeVelocity * theta[end] / nu)^(-0.268)

# Solving for the displacement thickness
displacementThickness = H.*theta

## Storing values to plot for Head's Method
theta_Head = theta
displacementThickness_Head = displacementThickness
H_Head = H
cf_Head = cf
Cf = sum(cf.*dx)
println("Head's Method Cf = ",Cf)

####################################################
# Schlichting Solution for turbulent boundary layers
Rex = edgeVelocity .* x ./ nu
ReL = edgeVelocity * plateLength / nu
displacementThickness_Schlicting = zeros(length(x))
theta_Schlicting = zeros(length(x))
cf_Schlicting = zeros(length(x))
for i = 1:length(x)
    displacementThickness_Schlicting[i] = 0.046 * x[i] / Rex[i]^(0.2)
    theta_Schlicting[i] = 0.036 * x[i] / Rex[i]^(0.2)
    cf_Schlicting[i] = 0.0592 / Rex[i]^(0.2)
end

Cf = 0.074 / ReL^(0.2)
println("Schlichting Solution Cf = ",Cf)

## Storing values to plot for Schlicting Solution
theta_Schlicting = theta_Schlicting
displacementThickness_Schlicting = displacementThickness_Schlicting
H_Schlicting = displacementThickness_Schlicting./theta_Schlicting
cf_Schlicting = cf_Schlicting

####################################################
# Plotting the various methods together

figure()
title("Theta")
xlabel("x (meters)")
ylabel("Height (meters)")
plot(x,theta_Thwaite,label = "Thwaite's Method",linestyle = "--",linewidth = 6)
plot(x,theta_Blasius,label = "Blasius Solution",linestyle = "-",linewidth = 3)
plot(x,theta_Head,label = "Head's Method",linestyle = "--",linewidth = 6)
plot(x,theta_Schlicting,label = "Schlichting Solution",linestyle = "-",linewidth = 3)
legend()
grid("on")

figure()
title("Displacement Thickness")
xlabel("x (meters)")
ylabel("Height (meters)")
plot(x,displacementThickness_Thwaite,label = "Thwaite's Method",linestyle = "--",linewidth = 6)
plot(x,displacementThickness_Blasius,label = "Blasius Solution",linestyle = "-",linewidth = 3)
plot(x,displacementThickness_Head,label = "Head's Method",linestyle = "--",linewidth = 6)
plot(x,displacementThickness_Schlicting,label = "Schlichting Solution",linestyle = "-",linewidth = 3)
legend()
grid("on")

figure()
title("Health Factor")
xlabel("x (meters)")
ylabel("Value")
plot(x,H_Thwaite,label = "Thwaite's Method",linestyle = "--",linewidth = 6)
plot(x,H_Blasius,label = "Blasius Solution",linestyle = "-",linewidth = 3)
plot(x,H_Head,label = "Head's Method",linestyle = "--",linewidth = 6)
plot(x,H_Schlicting,label = "Schlichting Solution",linestyle = "-",linewidth = 3)
legend()
grid("on")

figure()
title("Local Drag Coefficient")
xlabel("x (meters)")
ylabel("Value")
plot(x,cf_Thwaite,label = "Thwaite's Method",linestyle = "--",linewidth = 6)
plot(x,cf_Blasius,label = "Blasius Solution",linestyle = "-",linewidth = 3)
plot(x,cf_Head,label = "Head's Method",linestyle = "--",linewidth = 6)
plot(x,cf_Schlicting,label = "Schlichting Solution",linestyle = "-",linewidth = 3)
legend()
grid("on")
ylim(0,0.01)
