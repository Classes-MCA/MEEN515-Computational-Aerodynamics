# include("Box/Classes/Aerodynamics/Homework 4/Efficiency Convergence Plots.jl")

using Pkg
Pkg.add("PyPlot")
using PyPlot

numberOfPanels = [1,2,4,6,8,10,12,14,16,18,20,30,40,100,200,400,800]

#efficiency = [0.667,0.81,0.868,0.899,0.919,0.932,0.941,0.949,0.954,0.959,0.972,0.979,0.992,0.996,0.998,0.999]

efficiency = zeros(length(numberOfPanels))

for i = 1:length(numberOfPanels)

    println("numPanels = ",numberOfPanels[i])

    CL, CDi_far, CDi_near, efficiency[i] = VLM_old(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,alpha,numberOfPanels[i],showPanels,showFlow,useTips,tipHeight,showPlots)

    println("")

end

figure()
plot(numberOfPanels,efficiency)
scatter(numberOfPanels,efficiency)
title("Efficiency Convergence for Elliptic Lift Distribution")
xlabel("Number of Panels")
ylabel("Efficiency")
xscale("log")
#xlim(1,10^3)
grid("on")