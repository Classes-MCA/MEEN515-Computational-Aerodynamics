# This function takes a NACA airfoil and returns the inviscid solution
# for the flow using the Hess-Smith panel method.
# The pressure data from the Hess-Smith panel method is then used as
# the input for the edge of the boundary layer.

# Inputs:
#   e: max camber
#   p: max camber location
#   t: max thickness
#   dphi: angular grid spacing for cosine spacing
#   angle: angle of attack (degrees)
#   freestream: freestream velocity
#   density: free field density
#   showPlot: Boolean, whether to plot the pressure distribution

# include("Box/Classes/Aerodynamics/Homework 3/miniXFOIL.jl")

include("HessSmith.jl")

using Pkg
Pkg.add("PyPlot")
using PyPlot

# Mark, trust me. Take the time to plan this code out and it will go by sooooo much faster.
# Mark, just do it.

function flattenSurface(surfacePanels)

    controlPoints = zeros(length(surfacePanels[:,1]) - 1,2)
    zCoordinates  = zeros(length(surfacePanels[:,1]) - 2,1) # One less z-coordinate than control points

    for i = 1:length(controlPoints[:,1])

        xControlPoint = (surfacePanels[i,1] + surfacePanels[i+1,1])/2
        yControlPoint = (surfacePanels[i,2] + surfacePanels[i+1,2])/2

        controlPoints[i,:] = [xControlPoint,yControlPoint]

        #println(xControlPoint,' ',yControlPoint)

    end

    #figure()
    #scatter(controlPoints[:,1],controlPoints[:,2])

    for i = 2:length(zCoordinates)

        zCoordinates[i] = zCoordinates[i-1] + sqrt((controlPoints[i + 1,1] - controlPoints[i,1])^2 + (controlPoints[i + 1,2] - controlPoints[i,2])^2)

    end

    return zCoordinates

end

function interpolate(oldArray,newDensity)

    newArray = []

    for i = 1:length(oldArray)-1

        firstPoint  = oldArray[i]
        secondPoint = oldArray[i+1]

        difference = secondPoint - firstPoint

        for j = 0:newDensity-1

            # interpolatedPoint = (firstPoint + secondPoint) / newDensity
            # newArray[i+(i+j)] = interpolatedPoint

            interpolatedPoint = firstPoint + difference*j/newDensity
            append!(newArray,interpolatedPoint)
            
        end

    end

    # extrapolating final values
    for j = 0:newDensity-1

        firstPoint  = oldArray[end-1]
        secondPoint = oldArray[end]
        difference = secondPoint-firstPoint

        extrapolatedPoint = firstPoint + difference*j/newDensity

        append!(newArray,extrapolatedPoint)


    end

    #println(newArray)

    #println("length old array = ",length(oldArray))
    #println("length new array = ",length(newArray))

    return newArray

end

function determinePressureGradient(pressureCoefficient,zCoordinates,freestream,density)

    dynamicPressure = 0.5 * density * freestream^2
    pressure = pressureCoefficient.*dynamicPressure 
    #println(pressureCoefficient)
    #println(pressure)

    dPdz = zeros(length(zCoordinates)-1,1)

    # First point: forward difference
    dz = zCoordinates[2] - zCoordinates[1]
    dPdz[1] = (pressure[2] - pressure[1]) / dz

    # Rest of the points: centered difference
    for i = 2:length(dPdz)

        dz = zCoordinates[i + 1] - zCoordinates[i-1]
        dPdz[i] = (pressure[i+1] - pressure[i-1]) / dz

        # println("i = ",i," dz = ",dz," dPdz = ",dPdz[i])

    end

    return dPdz

end

function determineVelocityGrandient(dPdz,edgeVelocity,density,zCoordinates)

    dVedz = zeros(length(dPdz))
    
    # # First point: forward difference
    # dz = zCoordinates[2] - zCoordinates[1]
    # dVedz[1] = (edgeVelocity[2] - edgeVelocity[1])/dz

    # for i = 2:length(dVedz)

    #     #dVedz[i] = -dPdz[i]/edgeVelocity[i]/density # Old version

    #     dz = zCoordinates[i + 1] - zCoordinates[i-1]
    #     dVedz[i] = (edgeVelocity[i+1] - edgeVelocity[i-1]) / dz

    #     #dVedz[i]  = -dPdz[i]*edgeVelocity[i]/density/((edgeVelocity[i]/343)^2) # Speed of sound one
    #     #println(dVedz[i])
    #     #println("i = ",i," dVedz = ",dVedz[i])

    # end

    for i = 1:length(dVedz)

        dVedz[i] = -dPdz[i]/edgeVelocity[i]/density 

    end
    #println(dVedz)
    return dVedz

end

function LaminarBL(surface,edgeVelocity,dVedz,thetaInitial,nu)

    # Uses Thwaite's Method until laminar separation
    # I think this is working, tests with no pressure or velocity gradient work
    # So it's either a problem with the input or the process in here

    # Getting theta arrays ready
    theta = zeros(length(surface));
    theta[1] = thetaInitial;
    dThetadx = zeros(length(surface))

    # Other variables
    lambda = zeros(length(surface))
    l = zeros(length(surface))
    H = zeros(length(surface))
    displacementThickness = zeros(length(surface))
    cf = zeros(length(surface))
    transitionIndex = length(surface)

    for i = 1:length(surface) - 1

        # Getting our dz variable
        dz = surface[i + 1] - surface[i]
        # println(dz)

        # Calculate lambda
        lambda[i] = theta[i]^2 / nu * dVedz[i]
        #println(dVedz[i])

        # Solving the differential equation for this point
        dThetadx[i] = 0.225 * nu / (edgeVelocity[i] * theta[i]) - 3 * (theta[i] / edgeVelocity[i]) * dVedz[i]
        # println(dThetadx[i])
        # println(edgeVelocity[i])

        # Adding the solution to the next point
        theta[i + 1] = theta[i] + (dThetadx[i] * dz)
    
        # Calculating l and H
        if lambda[i] >= 0 && lambda[i] <= 0.1
            # println("0 <= ",lambda[i]," <= 0.1")
            l[i] = 0.22 + 1.57*lambda[i] - 1.8*lambda[i]^2
            H[i] = 2.61 - 3.75*lambda[i] - 5.24*lambda[i]^2
        elseif lambda[i] > -0.1 && lambda[i] < 0
            # println("-0.1 < ",lambda[i]," < 0")
            l[i] = 0.22 + 1.402*lambda[i] + 0.018*lambda[i] / (0.107 + lambda[i])
            H[i] = 2.088 + 0.0731 / (0.14 + lambda[i])
        elseif lambda[i] > 0.1
            # println(lambda[i], " > 0.1")
            l[i] = 0.22 + 1.57*0.1 - 1.8*0.1^2
            H[i] = 2.61 - 3.75*0.1 - 5.24*0.1^2
        elseif lambda[i] < -0.1
            println("Laminar Separation predicted at z = ", surface[i])
            return displacementThickness, theta, H, cf, transitionIndex
        end

        ## Calculate delta star (displacement thickness)
        displacementThickness[i] = H[i].*theta[i]

        ## Calculate cf (local skin friction coefficient)
        #cf[i] = 2 .* (dThetadx[i] + dVedz[i]*theta[i]/edgeVelocity[i] * (H[i] + 2))
        cf[i] = 2 * nu * l[i] / (theta[i] * edgeVelocity[i])

        # Checking for transition # FIXME: It only has transition if the freestream is sufficiently large, seems wrong
        if H[i] > 2.1 && H[i] < 2.8
            Rez = edgeVelocity[i] * surface[i] / nu
            checkRez =  -40.4557 + 64.8066*H[i] - 26.7538*H[i]^2 + 3.3819*H[i]^3
            #println("log10(Rez) = ",log10(Rez),"  checkRez = ",checkRez)
            if log10(Rez) > checkRez
                transitionIndex = i
                #println("transitioned at z = ",surface[i])
                println("-------------------------------------------")
                return displacementThickness, theta, H, cf, transitionIndex
            end
            
        end

        # Checking for transition another way (Michel's Method)
        Rez = edgeVelocity[i] * surface[i] / nu
        Retheta = edgeVelocity[i] * theta[i]   / nu

        if Retheta >= (1.174 * (Rez^0.46 + 22400*Rez^(-0.54)))
            println("New Method")
            transitionIndex = i
            return displacementThickness, theta, H, cf, transitionIndex
        end

    end

    return displacementThickness, theta, H, cf, transitionIndex

end

function TurbulentBL(surface,edgeVelocity,dVedz,thetaInitial,nu,transitionIndex)
    
    H = zeros(length(surface))
    H1 = zeros(length(surface))
    cf = zeros(length(surface))
    theta = zeros(length(surface))
    displacementThickness = zeros(length(surface))
    separationIndex = length(surface)

    H[transitionIndex] = 1.28 # From a flat plate, the notes said use this value
    theta[transitionIndex] = thetaInitial

    # Solving for the intial H1 value
    if H[transitionIndex] <= 1.6
        H1[transitionIndex] = 0.8234*(H[transitionIndex] - 1.1)^(-1.287) + 3.3
    else
        H1[transitionIndex] = 1.5501*(H[transitionIndex] - 0.6778)^(-3.064) + 3.3
    end

    for i = transitionIndex:length(surface)-1

        # Getting our dz variable
        dz = surface[i + 1] - surface[i]

        #println(H1[i])
        # Solving for H
        if H1[i] >= 5.3
            H[i] = 0.86*(H1[i] - 3.3)^(-0.777) + 1.1
        elseif H1[i] > 3.3
            H[i] = 1.1538*(H1[i] - 3.3)^(-0.326) + 0.6778
        else
            H[i] = 3.0
            #return displacementThickness, theta, H, cf, separationIndex
        end

        displacementThickness[i] = H[i]*theta[i]

        # Solving for cf
        cf[i] = 0.246*10^(-0.678*H[i]) * (edgeVelocity[i] * theta[i] / nu)^(-0.268)
    
        # Solving for dThetadx
        dThetadx = cf[i]/2 - dVedz[i]*theta[i]/edgeVelocity[i] * (H[i] + 2)
        
        # Solving for dH1dx
        if H1[i] >= 3
            dH1dx = 0.0306/theta[i] * (H1[i] - 3)^(-0.6169) - dVedz[i]*H1[i]/edgeVelocity[i] - dThetadx*H1[i]/theta[i]
        else
            dH1dx = 0
            #return displacementThickness, theta, H, cf, separationIndex
        end

        if H[i] >= 3.0
            separationIndex = i
            #println("Turbulent Separation at z = ",surface[separationIndex])
            return displacementThickness, theta, H, cf, separationIndex
        end

        # dz may become negative at the very end of the surface array
        theta[i + 1] = theta[i] + dThetadx*dz
        H1[i + 1]    = H1[i] + dH1dx*dz

    end

    # In case it never separates
    separationIndex = length(surface) - 1
    return displacementThickness, theta, H, cf, separationIndex

end

# TEST CASES
# Force the pressure gradient to be zero
# Force the velocity gradient to be zero
# Force the edge velocity to be constant

# BOTTOM SURFACE IS WORKING. DON'T YOU DARE TOUCH IT!!!

function miniXFOIL(e,p,t,dphi,angle,freestream,density,nu,showPlot)

    interpolationAmount = 1

    pressureCoefficient, CL_integration, CM_integration, x, panels, Velocity = HessSmith(e,p,t,dphi,angle,freestream,density,showPlot);
    #plot(x,pressureCoefficient[Int(length(pressureCoefficient)/2):end])
    # Need to find the stagnation point based on the angle of attack

    # Region to check for stagnation point include
    #lowerSurfaceStart = Int(floor(0.25*length(pressureCoefficient)))
    #upperSurfaceEnd   = Int(floor(0.75*length(pressureCoefficient)))
    stagnationPointIndex = findmax(pressureCoefficient)[2]


    # Defining the top and bottom surfaces and then adding a point on the plot that shows
    # the location of the stagnation point
    bottomSurfacePanels = reverse(panels[1:stagnationPointIndex,:],dims=1) # Reversed to show proper direction along chord
    topSurfacePanels    = panels[stagnationPointIndex:end,:]
    # println(bottomSurfacePanels[:,1])
    # println(topSurfacePanels[:,1])
    # plot(bottomSurfacePanels[:,1],bottomSurfacePanels[:,2])
    # plot(topSurfacePanels[:,1],topSurfacePanels[:,2])

    if showPlot == true 
        scatter(panels[stagnationPointIndex,1],panels[stagnationPointIndex,2],color = "red")
    end

    # Turn the top and bottom surfaces of the airfoil into flat plates by creating
    # a coordinate 'z' that lies on the surface of the airfoil
    bottomSurface = flattenSurface(bottomSurfacePanels) # Good
    topSurface    = flattenSurface(topSurfacePanels) # Good

    # Interpolating
    # KNOWN PROBLEM: The last two points in the new array will be duplicates of the second to last two points
    # bottomSurface = interpolate(bottomSurface,interpolationAmount)
    # topSurface    = interpolate(topSurface,interpolationAmount)
    # pressureCoefficient = interpolate(pressureCoefficient,interpolationAmount)
    # stagnationPointIndex = findmax(pressureCoefficient)[2] # Finding it again on the interpolated array
    # Velocity = interpolate(Velocity,interpolationAmount)
    #scatter(bottomSurface,zeros(length(bottomSurface)))
    #scatter(topSurface,ones(length(topSurface)))


    # determine the pressure gradient as a function of z
    dPdz_bottom = -determinePressureGradient(reverse(pressureCoefficient[1:stagnationPointIndex-1],dims=1),bottomSurface,freestream,density) # Forced to be positive at the beginning by using the "-"
    dPdz_top    = determinePressureGradient(pressureCoefficient[stagnationPointIndex+1:end],topSurface,freestream,density)

    # determine the velocity gradient as a function of z
    Ve_bottom = -reverse(Velocity[1:stagnationPointIndex-1],dims=1) # Made negative to give proper orientation. All but last element pulled because of a sign error
    Ve_top    = Velocity[stagnationPointIndex+1:end]

    dVedz_bottom = -determineVelocityGrandient(dPdz_bottom,Ve_bottom,density,bottomSurface)
    dVedz_top    = determineVelocityGrandient(dPdz_top,   Ve_top   ,density,topSurface)

    # Top surface
    # Laminar Section
    thetaInitial = sqrt(0.075*nu / dVedz_top[1])
    displacementThickness_top_Laminar, theta_top_Laminar, H_top_Laminar, cf_top_Laminar, transitionIndex_top = LaminarBL(topSurface,Ve_top,dVedz_top,thetaInitial,nu);

    # Turbulent Section
    thetaInitial = theta_top_Laminar[transitionIndex_top]
    displacementThickness_top_Turbulent, theta_top_Turbulent, H_top_Turbulent, cf_top_Turbulent, separationIndex_top = TurbulentBL(topSurface,Ve_top,dVedz_top,thetaInitial,nu,transitionIndex_top);
    
    # Bottom surface - FIXME: Not working
    # Laminar Section
    thetaInitial = sqrt(0.075*nu / dVedz_bottom[1])
    displacementThickness_bottom_Laminar, theta_bottom_Laminar, H_bottom_Laminar, cf_bottom_Laminar, transitionIndex_bottom = LaminarBL(bottomSurface,Ve_bottom,dVedz_bottom,thetaInitial,nu);

    # Turbulent Section
    thetaInitial = theta_bottom_Laminar[transitionIndex_bottom]
    displacementThickness_bottom_Turbulent, theta_bottom_Turbulent, H_bottom_Turbulent, cf_bottom_Turbulent, separationIndex_bottom = TurbulentBL(bottomSurface,Ve_bottom,dVedz_bottom,thetaInitial,nu,transitionIndex_bottom);
    
    ##########################################

    # Calculating Drag
    c = 1 # Chord Length
    drag_top = (2*theta_top_Turbulent[separationIndex_top] / c) * (Ve_top[separationIndex_top]/freestream)^((H_top_Turbulent[separationIndex_top] + 5)/2)
    drag_bottom = (2*theta_bottom_Turbulent[separationIndex_bottom] / c) * (Ve_bottom[separationIndex_bottom]/freestream)^((H_bottom_Turbulent[separationIndex_bottom] + 5)/2)
    Cd = drag_top + drag_bottom

    #Cd = drag / (0.5 * density * freestream^2 * c)
    println("Drag Coefficient = ",Cd)

    ##########################################

    # Displacement Thickness
    figure()
    plot(topSurface[1:transitionIndex_top],displacementThickness_top_Laminar[1:transitionIndex_top],label = "Laminar")
    plot(topSurface[transitionIndex_top:separationIndex_top],displacementThickness_top_Turbulent[transitionIndex_top:separationIndex_top],label = "Turbulent")
    scatter(topSurface[transitionIndex_top],displacementThickness_top_Laminar[transitionIndex_top],color = "red",label = "Transition")
    scatter(topSurface[separationIndex_top],displacementThickness_top_Turbulent[separationIndex_top],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Displacement Thickness - Top")
    legend()

    # theta
    figure()
    plot(topSurface[1:transitionIndex_top],theta_top_Laminar[1:transitionIndex_top],label = "Laminar")
    plot(topSurface[transitionIndex_top:separationIndex_top],theta_top_Turbulent[transitionIndex_top:separationIndex_top],label = "Turbulent")
    scatter(topSurface[transitionIndex_top],theta_top_Laminar[transitionIndex_top],color = "red",label = "Transition")
    scatter(topSurface[separationIndex_top],theta_top_Turbulent[separationIndex_top],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Theta - Top")
    legend()

    # H
    figure()
    plot(topSurface[1:transitionIndex_top],H_top_Laminar[1:transitionIndex_top],label = "Laminar")
    plot(topSurface[transitionIndex_top:separationIndex_top],H_top_Turbulent[transitionIndex_top:separationIndex_top],label = "Turbulent")
    scatter(topSurface[transitionIndex_top],H_top_Laminar[transitionIndex_top],color = "red",label = "Transition")
    scatter(topSurface[separationIndex_top],H_top_Turbulent[separationIndex_top],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("H - Top")
    legend()

    # cf
    figure()
    plot(topSurface[1:transitionIndex_top],cf_top_Laminar[1:transitionIndex_top],label = "Laminar")
    plot(topSurface[transitionIndex_top:separationIndex_top],cf_top_Turbulent[transitionIndex_top:separationIndex_top],label = "Turbulent")
    scatter(topSurface[transitionIndex_top],cf_top_Laminar[transitionIndex_top],color = "red",label = "Transition")
    scatter(topSurface[separationIndex_top],cf_top_Turbulent[separationIndex_top],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Skin Friction Coefficient - Top")
    legend()
    ylim(0,0.02)

    # # ------------------------

    # Displacement Thickness - I HACKED THE SEPARATION INDEX TO PUT IT WHERE I WANTED IT. SEE END OF TURBULENT BL FUNCTION
    figure()
    plot(bottomSurface[1:transitionIndex_bottom],displacementThickness_bottom_Laminar[1:transitionIndex_bottom],label = "Laminar")
    plot(bottomSurface[transitionIndex_bottom:separationIndex_bottom],displacementThickness_bottom_Turbulent[transitionIndex_bottom:separationIndex_bottom],label = "Turbulent")
    scatter(bottomSurface[transitionIndex_bottom],displacementThickness_bottom_Laminar[transitionIndex_bottom],color = "red",label = "Transition")
    scatter(bottomSurface[separationIndex_bottom],displacementThickness_bottom_Turbulent[separationIndex_bottom],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Displacement Thickness - Bottom")
    legend()

    # theta
    figure()
    plot(bottomSurface[1:transitionIndex_bottom],theta_bottom_Laminar[1:transitionIndex_bottom],label = "Laminar")
    plot(bottomSurface[transitionIndex_bottom:separationIndex_bottom],theta_bottom_Turbulent[transitionIndex_bottom:separationIndex_bottom],label = "Turbulent")
    scatter(bottomSurface[transitionIndex_bottom],theta_bottom_Laminar[transitionIndex_bottom],color = "red",label = "Transition")
    scatter(bottomSurface[separationIndex_bottom],theta_bottom_Turbulent[separationIndex_bottom],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Theta - Bottom")
    legend()

    # H
    figure()
    plot(bottomSurface[1:transitionIndex_bottom],H_bottom_Laminar[1:transitionIndex_bottom],label = "Laminar")
    plot(bottomSurface[transitionIndex_bottom:separationIndex_bottom],H_bottom_Turbulent[transitionIndex_bottom:separationIndex_bottom],label = "Turbulent")
    scatter(bottomSurface[transitionIndex_bottom],H_bottom_Laminar[transitionIndex_bottom],color = "red",label = "Transition")
    scatter(bottomSurface[separationIndex_bottom],H_bottom_Turbulent[separationIndex_bottom],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("H - Bottom")
    legend()

    # cf
    figure()
    plot(bottomSurface[1:transitionIndex_bottom],cf_bottom_Laminar[1:transitionIndex_bottom],label = "Laminar")
    plot(bottomSurface[transitionIndex_bottom:separationIndex_bottom],cf_bottom_Turbulent[transitionIndex_bottom:separationIndex_bottom],label = "Turbulent")
    scatter(bottomSurface[transitionIndex_bottom],cf_bottom_Laminar[transitionIndex_bottom],color = "red",label = "Transition")
    scatter(bottomSurface[separationIndex_bottom],cf_bottom_Turbulent[separationIndex_bottom],color = "purple",label = "Separation")
    xlabel("Distance Along Airfoil Surface")
    ylabel("Magnitude")
    grid("on")
    title("Skin Friction Coefficient - Bottom")
    legend()
    ylim(0,0.02)

end