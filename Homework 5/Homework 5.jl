# Purpose: model a supersonic airfoil using supersonic thin airfoil theory and Shock-Expansion theory

# include("Box/Classes/Aerodynamics/Homework 5/Homework 5.jl")

using Pkg
Pkg.add("PyPlot")
using PyPlot

# Define univeral variables
AOA = -5:10
alpha = AOA.*pi/180 # converting to radians
machNumber = 1.4
chord = 1
maxThickness = 0.05*chord
averageCamber = 0
gamma = 1.4

## Thin Airfoil theory

# Lift Coefficient
cl_taf = 4 .* alpha./sqrt(machNumber^2 - 1)

# Drag Coefficient
cd_taf = 4/sqrt(machNumber^2 - 1) .* (alpha.^2 .+ (maxThickness/chord)^2)

# Pitching Moment about the leading edge
cm_taf = -2/sqrt(machNumber^2 - 1) .* (alpha .+ averageCamber/chord)


## Shock-Expansion theory

function shock(theta,M_1)

    # Solve numerically for Beta
    Beta = 1 * pi/180 # initial guess
    iter = 0
    maxIter = 1000
    error = 1
    maxError = 10^(-6)
    weightFactor = 3 # How much weight to put on the rootFunction in the solver. Speeds up solver but risks overshooting answer
    while abs(error) > maxError

        # eqn 9.23
        LHS = tan(theta)
        RHS = 2*cot(Beta) * (M_1^2*sin(Beta)^2 - 1) / (M_1^2 * (gamma + cos(2*Beta)) + 2)
        rootFunction = LHS - RHS

        if rootFunction > maxError
            #println("Beta on iteration ",iter," = ",Beta*180/pi)
            Beta = Beta + 0.5 * pi/180 * rootFunction * weightFactor
        else
            #println("Breaking...")
            break
        end

        iter = iter + 1
        if iter > maxIter
            println("Beta did not converge in shock() function")
            break
        end

    end

    # Solving for M_2
    Mn_1 = M_1*sin(Beta) # eqn 9.13
    Mn_2 = sqrt((1 + 0.5*(gamma - 1)*Mn_1^2) / (gamma*Mn_1^2 - 0.5*(gamma-1))) # eqn 9.14
    #println("Beta = ",Beta*180/pi)
    M_2 = Mn_2 / sin(Beta - theta) # eqn 9.18

    # Solving for the pressure pressure ratio
    pressureRatio = 1 + (2*gamma)/(gamma+1) * (Mn_1^2 - 1) # eqn 9.16

    # println("Theta = ",theta*180/pi)
    # println("Beta = ",Beta*180/pi)
    # println("Mn_1 = ",Mn_1)
    # println("Mn_2 = ",Mn_2)
    # println("M_2 = ",M_2)
    # println("P2/P1 = ",pressureRatio)

    return pressureRatio,M_2

end

function expansion(theta,M_1)

    nu_M_1 = prandtlMeyer(gamma,M_1) # eqn 9.42
    nu_M_2 = theta + nu_M_1
    # println(nu_M_2*180/pi)
    M_2 = reversePrandtlMeyer(gamma,nu_M_2)

    pressureRatio = ((1 + 0.5*(gamma - 1)*M_1^2) / (1 + 0.5*(gamma - 1)*M_2^2))^(gamma/(gamma - 1))

    return pressureRatio,M_2
end

function prandtlMeyer(gamma,M)

    nu = sqrt((gamma + 1) / (gamma - 1))
    nu = nu * atan(sqrt((gamma - 1)/(gamma + 1) * (M^2 - 1)))
    nu = nu - atan(sqrt(M^2 - 1))
    #nu = nu * 180/pi # Convert to degrees

    return nu
end

function reversePrandtlMeyer(gamma,nu_M)

    # Solve numerically for M
    M = 1.2 # initial guess
    iter = 0
    maxIter = 1000
    error = 1
    maxError = 10^(-6)
    weightFactor = 1
    while abs(error) > maxError

        # eqn 9.42
        LHS = nu_M
        RHS = prandtlMeyer(gamma,M)
        rootFunction = LHS - RHS

        if rootFunction > maxError
            #println("M on iteration ",iter," = ",M)
            M = M + 0.05 * rootFunction * weightFactor
        else
            #println("Breaking...")
            break
        end

        iter = iter + 1
        if iter > maxIter
            println("M did not converge in reversePrandtlMeyer() function")
            break
        end

    end

    return M
end

function defineDiamondAirfoil(chord,maxThickness)

    airfoilCoordinates = zeros(4,2) # Initializing the coordiantes array

    # Leading edge
    airfoilCoordinates[1,1] = 0
    airfoilCoordinates[1,2] = 0

    # Center of upper surface
    airfoilCoordinates[2,1] = chord/2
    airfoilCoordinates[2,2] = maxThickness/2 * chord

    # Trailing edge
    airfoilCoordinates[3,1] = chord
    airfoilCoordinates[3,2] = 0

    # Center of lower surface
    airfoilCoordinates[4,1] = chord/2
    airfoilCoordinates[4,2] = -maxThickness/2 * chord

    return airfoilCoordinates
end

#-- Main code block

# Initializing Arrays
CL = zeros(length(alpha))
CD = zeros(length(alpha))
CM = zeros(length(alpha))

# Define the airfoil
airfoilCoordinates = defineDiamondAirfoil(chord,maxThickness)

# Panel numbering is: (Dr. Ning's cm equation will need to be adapted)
# 1 = leading panel on upper surface
# 2 = trailing panel on upper surface
# 3 = leading panel on lower surface
# 4 = trailing panel on lower surface

# Looping through all the angles of attack
for i = 1:length(alpha)

    println("Alpha = ",alpha[i]*180/pi)

    # Panel 1
    # Determine whether the angle of attack is larger than the panel angle
    geometricAngle = atan(airfoilCoordinates[2,2],airfoilCoordinates[2,1])
    # println("Geometric Angle = ",geometricAngle*180/pi)
    theta_1 = -(alpha[i] - geometricAngle)
    if theta_1 > 0
        pressureRatio_P1_Pinf,M_1 = shock(theta_1,machNumber)
    else
        theta_1 = abs(theta_1) # preserving positive value
        pressureRatio_P1_Pinf,M_1 = expansion(theta_1,machNumber)
    end

    # Panel 2
    # Assume expansion wave
    theta_2 = 2*atan(airfoilCoordinates[2,2],airfoilCoordinates[2,1])
    pressureRatio_P2_P1,M_2 = expansion(theta_2,M_1)

    # Panel 3
    # Determine whether the angle of attack is larger than the panel angle
    theta_3 = alpha[i] + geometricAngle
    if theta_3 > 0
        pressureRatio_P3_Pinf,M_3 = shock(theta_3,machNumber)
    else
        theta_3 = abs(theta_3) # preserving positive value
        pressureRatio_P3_Pinf,M_3 = expansion(theta_3,machNumber)
    end

    # Panel 4
    # Assume expansion wave
    theta_4 = 2*atan(airfoilCoordinates[2,2],airfoilCoordinates[2,1]) # Symmetric with Panel 2
    pressureRatio_P4_P3,M_4 = expansion(theta_4,M_3)

    # Now we will use the pressure distributions to calculate lift, drag, and moment coefficients
    # See Example 9.12 in Anderson's book
    # Lift = Pressure Difference
    # CL = P_bottom - P_top / (q*S)
    lift_1 = -pressureRatio_P1_Pinf * cos(geometricAngle - alpha[i]) * cos(alpha[i]) * chord/2
    lift_2 = -pressureRatio_P2_P1 * pressureRatio_P1_Pinf * cos(-geometricAngle - alpha[i]) * cos(alpha[i]) * chord/2
    lift_3 = pressureRatio_P3_Pinf * cos(geometricAngle + alpha[i]) * cos(alpha[i]) * chord/2
    lift_4 = pressureRatio_P4_P3 * pressureRatio_P3_Pinf * cos(-geometricAngle + alpha[i]) * cos(alpha[i]) * chord/2
    CL[i] = (lift_1 + lift_2 + lift_3 + lift_4) / (0.5*gamma * machNumber^2 * chord)

    drag_1 = pressureRatio_P1_Pinf * sin(geometricAngle - alpha[i]) * cos(alpha[i]) * chord/2
    drag_2 = pressureRatio_P2_P1 * pressureRatio_P1_Pinf * sin(-geometricAngle - alpha[i]) * cos(alpha[i]) * chord/2
    drag_3 = pressureRatio_P3_Pinf * sin(geometricAngle + alpha[i]) * cos(alpha[i]) * chord/2
    drag_4 = pressureRatio_P4_P3 * pressureRatio_P3_Pinf * sin(-geometricAngle + alpha[i]) * cos(alpha[i]) * chord/2
    CD[i] = (drag_1 + drag_2 + drag_3 + drag_4) / (0.5*gamma * machNumber^2 * chord)

    moment_1 = pressureRatio_P1_Pinf / (4*gamma*machNumber^2) * (1 + (maxThickness/chord)^2)
    moment_2 = pressureRatio_P2_P1 * pressureRatio_P1_Pinf / (4*gamma*machNumber^2) * (3 - (maxThickness/chord)^2)
    moment_3 = pressureRatio_P3_Pinf / (4*gamma*machNumber^2) * (-1 - (maxThickness/chord)^2)
    moment_4 = pressureRatio_P4_P3 * pressureRatio_P3_Pinf / (4*gamma*machNumber^2) * (-3 + (maxThickness/chord)^2)
    CM[i] = moment_1 + moment_2 + moment_3 + moment_4

    # println("Angle of attack = ",AOA[i])

    # println("Lift_1 = ",lift_1)
    # println("Lift_2 = ",lift_2)
    # println("Lift_3 = ",lift_3)
    # println("Lift_4 = ",lift_4)

    # println(" ")

    # println("Drag_1 = ",drag_1)
    # println("Drag_2 = ",drag_2)
    # println("Drag_3 = ",drag_3)
    # println("Drag_4 = ",drag_4)

    # println("")

    # println("Moment_1 = ",moment_1)
    # println("Moment_2 = ",moment_2)
    # println("Moment_3 = ",moment_3)
    # println("Moment_4 = ",moment_4)

end

figure()
plot(AOA,cl_taf,label = "Thin Airfoil Theory")
scatter(AOA,CL,label = "Shock-Expansion Theory",color = "red")
title("Lift Coefficient vs AOA")
xlabel("Angle of Attack (Degrees)")
ylabel("Lift Coefficient")
grid("on")
legend()

figure()
plot(AOA,cd_taf,label = "Thin Airfoil Theory")
scatter(AOA,CD,label = "Shock-Expansion Theory",color = "red")
title("Drag Coefficient vs AOA")
xlabel("Angle of Attack (Degrees)")
ylabel("Drag Coefficient")
grid("on")
legend()

figure()
plot(AOA,cm_taf,label = "Thin Airfoil Theory")
scatter(AOA,CM,label = "Shock-Expansion Theory",color = "red")
title("Moment Coefficient about Leading Edge vs AOA")
xlabel("Angle of Attack (Degrees)")
ylabel("Moment Coefficient")
grid("on")
legend()