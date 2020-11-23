# Purpose: To implement a simple Blade Element Momentum theory code

# include("Box/Classes/Aerodynamics/Homework 6/Homework 6.jl")


using Pkg
Pkg.add("PyPlot")
Pkg.add("CSV")
Pkg.add("Interpolations")
Pkg.add("FLOWMath")
using PyPlot
using CSV
using Interpolations
using FLOWMath

CL_file = "Box/Classes/Aerodynamics/Homework 6/Airfoil Data CL.csv"
CD_file = "Box/Classes/Aerodynamics/Homework 6/Airfoil Data CD.csv"
Geometry_file = "Box/Classes/Aerodynamics/Homework 6/Propeller Geometry.csv"

RPMvals = [4005,6005,6710]

colors = ["red" "blue";"orange" "purple";"black" "green"]

for l = 1:length(RPMvals)

    RPM = RPMvals[l]

    if RPM == 4005
        Comparison_file = "Box/Classes/Aerodynamics/Homework 6/Comparison Data 4005 RPM.csv"
    elseif RPM == 6005
        Comparison_file = "Box/Classes/Aerodynamics/Homework 6/Comparison Data 6005 RPM.csv"
    elseif RPM == 6710
        Comparison_file = "Box/Classes/Aerodynamics/Homework 6/Comparison Data 6710 RPM.csv"
    else
        println("No comparison file matching that RPM value")
    end

    # Create a function that returns CL for a given angle of attack and mach number
    function getCL(alpha,machNumber,CL_file)

        CLdata = CSV.read(CL_file) # Downloading data
        angles = CLdata[:,1]       # Getting the angles array
        liftCoefficients =CLdata[:,2] # Getting the lift coefficients array
        cl = linear(angles,liftCoefficients,alpha) # Doing a linear interpolation (alpha in degrees)

        cl = cl / sqrt(1 - machNumber^2)

        return cl

    end

    # Create a function that returns CD for a given angle of attack and mach number
    function getCD(alpha,machNumber,CD_file)

        CDdata = CSV.read(CD_file) # Downloading data
        angles = CDdata[:,1]       # Getting the angles array
        dragCoefficients =CDdata[:,2] # Getting the drag coefficients array
        cd = linear(angles,dragCoefficients,alpha) # Doing a linear interpolation (alpha in degrees)

        #cd = cd / sqrt(1 - machNumber^2)

        return cd

    end

    # Create a residual function to calculate the residual of phi
    function Residual_phi(phi,theta,radialLocation,totalRadius,chordLength,numBlades,Omega,freestream,returnVals)

        machNumber = sqrt((Omega*radialLocation)^2 + freestream^2)/343 # speed of sound = 343 m/s

        # Using Radians
        theta = theta * pi/180
        alpha = theta - phi # Calculating the effective angle of attack

        # Getting the lift and drag coefficients (convert alpha to degrees)
        cl = getCL(alpha*180/pi,machNumber,CL_file)
        cd = getCD(alpha*180/pi,machNumber,CD_file)

        # Calculating the normal and tangential force components
        cn = cl * cos(phi) - cd * sin(phi)
        ct = cl * sin(phi) + cd * cos(phi)

        # Calculating the loss factor in section 6.1.2 in Dr. Ning's book
        Rhub = 0.010*totalRadius
        ftip = numBlades / 2 * (totalRadius - radialLocation) / (radialLocation * abs(sin(phi)))
        Ftip = 2/pi * acos(exp(-ftip))
        fhub = numBlades / 2 * (radialLocation - Rhub) / (Rhub*abs(sin(phi)))
        Fhub = 2/pi * acos(exp(-fhub))
        F = Ftip*Fhub # loss factor
        #println("r/R = ",radialLocation/totalRadius,". F = ",F)

        # Calculating the nondimensional factor sigmaPrime
        sigmaPrime = numBlades * chordLength / (2 * pi * radialLocation)

        # Calculating the 'a' value
        a = sigmaPrime * cn / (4 * F * sin(phi)^2 - sigmaPrime * cn)

        # Calculating tje 'a prime' value
        aPrime = sigmaPrime * ct / (4 * F * sin(phi) * cos(phi) + sigmaPrime * ct)

        if returnVals

            # Adjusting the normal force coefficient for the rotational correction, see eqn 9 in the paper Dr. Ning shared
            if radialLocation / totalRadius <= 0.85
                Vl = sqrt((Omega*radialLocation)^2 + freestream^2)
                cn = cn + 1.5 * (chordLength / radialLocation)^2 * (2*pi*alpha - cl) * (Omega * radialLocation / Vl)^2
            end

            return cn,ct,sigmaPrime,a,aPrime
        end

        return sin(phi) / (1 + a) - freestream / (Omega * radialLocation) * cos(phi) / (1 - aPrime)

    end

    # # Creating the plots for the first part of the Homework. It works!
    # angleRange = -30:0.5:30
    # liftCurve = zeros(length(angleRange))
    # dragCurve = zeros(length(angleRange))
    # machNumber = 0
    # for i = 1:length(angleRange)

    #     liftCurve[i] = getCL(angleRange[i],machNumber,CL_file)
    #     dragCurve[i] = getCD(angleRange[i],machNumber,CD_file)

    # end

    # # Creating the plots
    # figure()
    # plot(angleRange,liftCurve)
    # title("Lift Curve")
    # xlabel("Angle of Attack")
    # ylabel("CL")
    # grid("on")

    # figure()
    # semilogy(angleRange,dragCurve)
    # title("Drag Curve")
    # xlabel("Angle of Attack")
    # ylabel("CD")
    # grid("on")

    # Doing a loop over different advance ratios

    # Pulling in the comparison data
    comparisonData = CSV.read(Comparison_file)
    advanceRatios = comparisonData[:,1]
    CT_comparison = comparisonData[:,2]
    CP_comparison = comparisonData[:,3]
    CQ_comparison = zeros(length(CP_comparison))
    efficiency_comparison = comparisonData[:,4]

    CT = zeros(length(advanceRatios))
    CQ = similar(CT)
    CP = similar(CT)
    efficiency = similar(CT)

    for k = 1:length(advanceRatios[1:end])
        # Read in the geometry
        Radius = 0.127 # 5 inches
        geometryData = CSV.read(Geometry_file) # Downloading data
        givenRadialLocations = geometryData[:,1] * Radius # Radial Locations
        givenChordLengths = geometryData[:,2] * Radius # Chord Lengths
        givenTwistAngles = geometryData[:,3] # Twist angles

        radialLocations = givenRadialLocations[1]:(Radius-givenRadialLocations[1])/50:Radius
        chordLengths = linear(givenRadialLocations,givenChordLengths,radialLocations)
        twistAngles  = linear(givenRadialLocations,givenTwistAngles,radialLocations)

        # Define other parameters
        Omega = 2*pi*RPM/60 # eqn 6.65
        diameter = 2 * Radius
        freestream = advanceRatios[k] * RPM/60 * diameter # eqn 6.67
        rho = 1.225 # kg/m^3
        numBlades = 2

        # Initializing Arrays
        Thrust = zeros(length(radialLocations))
        Torque = similar(Thrust)
        Power  = similar(Thrust)

        for i = 1:length(radialLocations)

            #println("Advance Ratio: ",advanceRatios[k],". r/R = ",radialLocations[i]/Radius)

            lowerPhi = 0.01
            upperPhi = pi/2
            lowerBound = Residual_phi(lowerPhi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,false)
            upperBound = Residual_phi(upperPhi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,false)
            iter = 0
            residual = 1
            phi = 1 # initializing to avoid pass-by-reference
            while abs(residual) > 10^-6

                phi = (lowerPhi + upperPhi) / 2

                residual = Residual_phi(phi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,false)

                if sign(residual) == sign(lowerBound)
                    lowerPhi = phi
                    lowerBound = Residual_phi(lowerPhi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,false)
                else
                    upperPhi = phi
                    upperBound = Residual_phi(upperPhi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,false)
                end

                iter = iter + 1
                if iter > 100
                    println("Failed to converge")
                    break
                end

            end

            # println("phi = ",phi * 180/pi)

            # Calculating all the stuff to go with iter

            cn,ct,sigmaPrime,a,aPrime = Residual_phi(phi,twistAngles[i],radialLocations[i],Radius,chordLengths[i],numBlades,Omega,freestream,true)

            # Inflow Velocity (W) (convert twisting angle to radians)
            inflowVelocitySquared = (freestream * (1 + a))^2 + (Omega * radialLocations[i] * (1 - aPrime))^2
            #inflowVelocitySquared = (freestream * (1 + a) / sin(phi))^2
            #inflowVelocitySquared = (Omega * radialLocations[i] * (1 - aPrime)/cos(phi))^2

            # Thrust
            Thrust[i] = numBlades * cn * 0.5 * rho * inflowVelocitySquared * chordLengths[i]
            
            # Torque
            Torque[i] = numBlades * radialLocations[i] * ct * 0.5 * rho * inflowVelocitySquared * chordLengths[i]

        end

        # Creating coefficients
        totalThrust = trapz(radialLocations,Thrust)
        totalTorque = trapz(radialLocations,Torque)
        power = totalTorque * Omega

        n = Omega / (2*pi)
        D = diameter
        CT[k] = totalThrust / (rho * n^2 * D^4)
        CQ[k] = totalTorque / (rho * n^2 * D^5)
        CP[k] = power / (rho * n^3 * D^5)
        efficiency[k] = advanceRatios[k] * CT[k] / CP[k]

        CQ_comparison[k] = CP_comparison[k] ./ Omega .* n

        # # Plotting the current thrust distribution
        # figure()
        # plot(radialLocations/Radius,Thrust)
        # title(string("Thrust for Advance Ratio ",advanceRatios[k]))
        # xlabel("r/R")
        # ylabel("Thrust")

        # # Plotting the current torque distribution
        # figure()
        # plot(radialLocations/Radius,Torque)
        # title(string("Torque for Advance Ratio ",advanceRatios[k]))
        # xlabel("r/R")
        # ylabel("Torque")

    end

    # Plotting the final results on separate plots
    figure()
    scatter(advanceRatios,CT_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,CT,label = string("Calculated: ",RPM),color = colors[l,2])
    title(string("RPM = ",RPM))
    xlabel("Advance Ratio")
    ylabel("CT")
    ylim(-0.02,0.1)
    legend(loc = 3)
    grid("on")

    # figure()
    # scatter(advanceRatios,CQ_comparison,label = "UIUC")
    # scatter(advanceRatios,CQ,label = "Calculated")
    # title(string("RPM = ",RPM))
    # xlabel("Advance Ratio")
    # ylabel("CQ")
    # legend(loc = 3)
    # grid("on")

    figure()
    scatter(advanceRatios,CP_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,CP,label = string("Calculated: ",RPM),color = colors[l,2])
    title(string("RPM = ",RPM))
    xlabel("Advance Ratio")
    ylabel("CP")
    ylim(0,0.05)
    legend(loc = 3)
    grid("on")

    figure()
    scatter(advanceRatios,efficiency_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,efficiency,label = string("Calculated: ",RPM),color = colors[l,2])
    title(string("RPM = ",RPM))
    xlabel("Advance Ratio")
    ylabel("Efficiency")
    ylim(0,1)
    legend(loc = 3)
    grid("on")

    # Plotting the final results on the same plots
    figure(4)
    scatter(advanceRatios,CT_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,CT,label = string("Calculated: ",RPM),color = colors[l,2])
    title("Thrust Coefficient")
    xlabel("Advance Ratio")
    ylabel("CT")
    ylim(-0.02,0.1)
    legend(loc = 3)
    grid("on")

    figure(5)
    scatter(advanceRatios,CP_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,CP,label = string("Calculated: ",RPM),color = colors[l,2])
    title("Power Coefficient")
    xlabel("Advance Ratio")
    ylabel("CP")
    ylim(-0.02,0.05)
    legend(loc = 3)
    grid("on")

    figure(6)
    scatter(advanceRatios,efficiency_comparison,label = string("UIUC: ",RPM),color = colors[l,1])
    plot(advanceRatios,efficiency,label = string("Calculated: ",RPM),color = colors[l,2])
    title("Efficiency")
    xlabel("Advance Ratio")
    ylabel("Efficiency")
    ylim(0,1)
    legend()
    grid("on")

end