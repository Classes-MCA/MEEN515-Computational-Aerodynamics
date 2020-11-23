# This code is my panel method code

# include("Desktop/School/Aerodynamics/Homework 2/Homework 2-2.jl")

using Pkg
Pkg.add("PyPlot")
Pkg.add("LinearAlgebra")
using PyPlot
using LinearAlgebra

include("taf.jl") # We just need the naca4 function

# Function (Use): Make an x array that has cosine spacing
function makeXarray(c,dphi)

    phi = 0:dphi:pi
    x = zeros(length(phi),1)
    for i = 1:length(phi)
        x[i] = c/2 * (1-cos(phi[i]))
    end

    return x

end

# Function (Use): Create NACA Airfoil
function createNACA(e,p,t,x)

    T, ybar, dTdx, dybardx = naca4(e,p,t,x)

    return T, ybar

end

# Function (Use): Plot the airfoil
function plotAirfoil(panels,withTips)

    # Plot the airfoil
    figure()
    if withTips == true
        scatter(panels[:,1],panels[:,2],label = "Panel Tips")
    end
    plot(panels[:,1],panels[:,2],label = "Panel Surface")
    ylim(-.5,.5)
    grid("on")
    ylabel("Y")
    xlabel("X")
    title("Airfoil")
    legend()

end

# Function (Use): Plot the airfoil with control Points
function plotControlPoints(panels,controlPoints)

    plotAirfoil(panels)
    scatter(controlPoints[:,1],controlPoints[:,2],color = "red",label = "Control Points")
    #legend()

end

# Function (Use): Create the Panels
function createPanels(x,ybar,Thickness)

    panels = zeros(length(x)*2 - 1,2) # "-1" is so we can later get rid of the redunant (0,0) point

    for i = 1:(length(x)*2 - 1)

        if i <= length(x) # Fill in the bottom of the airfoil. Notice x is reversed to make panel always on right side
            panels[i,1] = x[length(x) - i + 1]
            panels[i,2] = ybar[length(x) - i + 1] - 0.5*Thickness[length(x) - i + 1]
        else # Fill in the top of the airfoil
            # The " + 1 " parts of this avoid having two points at the origin, which would result in a control point there
            panels[i,1] = x[i-length(x) + 1]
            panels[i,2] = ybar[i-length(x) + 1] + 0.5*Thickness[i-length(x) + 1]
        end

        #println(panels[i,1],' ',panels[i,2]) # Make sure panels are in proper order

    end

    return panels

end

# Function (Use): Calculate Panel Sines and Cosines of Angle
function calculateSineAndCosineOld(panels,numPanels)

    # Reference textbook by Dr. Ning for these equations

    panelSines = zeros(numPanels)
    panelCosines = zeros(numPanels)
    panelLengths = zeros(numPanels)

    for i = 1:numPanels

        deltaX = panels[i + 1,1] - panels[i,1]
        deltaY = panels[i + 1,2] - panels[i,2]
        length = sqrt(deltaX^2 + deltaY^2)
        panelLengths[i] = length
        panelSines[i] = deltaY/length
        panelCosines[i] = deltaX/length

    end

    return panelSines, panelCosines, panelLengths

end

function calculateSineAndCosine(currentPanel)

    # currentPanel = [X1 Y1; X2 Y2]

    deltaX = currentPanel[2][1,1] - currentPanel[1][1,1]
    #println(deltaX)
    deltaY = currentPanel[2][2,1] - currentPanel[1][2,1]
    #println(deltaY)
    length = sqrt(deltaX^2 + deltaY^2)
    panelLength = length
    panelSine = deltaY/length
    panelCosine = deltaX/length

    #println(currentPanel)
    #println(panelSine)

    return panelSine, panelCosine, panelLength

end

# Function (Use): Calculate panel control point location
function getControlPoint(firstCoordinate,secondCoordiante)

    #println(firstCoordinate,' ', secondCoordiante)

    # firstCoordinate = first coordiante of the panel (j)
    # secondCoordiante = second coordinate of the panel (j + 1), keeping airfoil on the right

    X_controlPoint = (firstCoordinate[1] + secondCoordiante[1])/2
    Y_controlPoint = (firstCoordinate[2] + secondCoordiante[2])/2

    # println(X_controlPoint,' ',Y_controlPoint)

    return [X_controlPoint,Y_controlPoint]

end

# Function (Use): Calculate unit normal to panel
function calculateUnitNormals(panels,numPanels,panelSines,panelCosines)

    unitNormals = zeros(numPanels,2)
    for i = 1:numPanels

        unitNormals[i,1] = -panelSines[i]
        unitNormals[i,2] = panelCosines[i]
        # println(sqrt(unitNormals[i,1]^2 + unitNormals[i,2]^2)) # Make sure they are unit vectors

    end

    return unitNormals
    
end

# Function (Use): Calculate unit tangent to panel
function calculateUnitTangents(panels,numPanels,panelSines,panelCosines)

    unitTangents = zeros(numPanels,2)
    for i = 1:numPanels

        # FIXME: Double-check this math
        unitTangents[i,1] = panelCosines[i]
        unitTangents[i,2] = panelSines[i]

        # println(sqrt(unitNormals[i,1]^2 + unitNormals[i,2]^2)) # Make sure they are unit vectors

    end

    return unitTangents
    
end

# Function (Use): Plot unit normal vectors
function plotUnitNormals(controlPoints,unitNormals,numPanels,scaleFactor)

    for i = 1:numPanels
        X = [controlPoints[i,1],(controlPoints[i,1] + unitNormals[i,1]./scaleFactor)]
        Y = [controlPoints[i,2],(controlPoints[i,2] + unitNormals[i,2]./scaleFactor)]

        plot(X,Y,color = "green")

    end

end

# Function (Use): Plot unit tangent vectors
function plotUnitTangents(controlPoints,unitTangents,numPanels,scaleFactor)

    for i = 1:numPanels
        X = [controlPoints[i,1],(controlPoints[i,1] + unitTangents[i,1]./scaleFactor)]
        Y = [controlPoints[i,2],(controlPoints[i,2] + unitTangents[i,2]./scaleFactor)]

        plot(X,Y,color = "purple")

    end

end


# Function (Use): Compute distances from both ends of a panel to a control point on another panel
function calculateDistances(currentPanel,currentControlPoint)

    # currentPanel = [X1 Y1; X2 Y2]
    # currentControlPoint = [Xc Yc]
    deltaX = currentPanel[1][1,1] - currentControlPoint[1]
    deltaY = currentPanel[1][2,1] - currentControlPoint[2]
    rij = sqrt(deltaX^2 + deltaY^2)

    #plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    deltaX = currentPanel[2][1,1] - currentControlPoint[1]
    deltaY = currentPanel[2][2,1] - currentControlPoint[2]
    rijPlusOne = sqrt(deltaX^2 + deltaY^2)

    #plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    return rij, rijPlusOne

end

# Function (Use): Make sure that the distances calculated are correct
function plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    X = [currentControlPoint[1],currentControlPoint[1] + deltaX]
    Y = [currentControlPoint[2],currentControlPoint[2] + deltaY]

    plot(X,Y,color = "black")

end

# Function (Use): Compute the angle Beta
function calculateBeta(currentPanel,currentControlPoint,i,j)

    # currentPanel = [X1 Y1; X2 Y2]
    # currentControlPoint = [Xc Yc]

    xj = currentPanel[1][1,1]
    xjPlus1 = currentPanel[2][1,1]

    #println("x",j," = ",xj,". x",j+1,"= ",xjPlus1)

    yj = currentPanel[1][2,1]
    yjPlus1 = currentPanel[2][2,1]
    #println("y",j," = ",yj,". y",j+1,"= ",yjPlus1)

    xi = currentControlPoint[1]
    yi = currentControlPoint[2]

    if i != j
        upperArgument = (xj - xi) * (yjPlus1 - yi) - (yj - yi) * (xjPlus1 - xi)
        lowerArgument = (xj - xi) * (xjPlus1 - xi) + (yj - yi) * (yjPlus1 - yi) # Should be '+' from dot product?
        Beta = atan(upperArgument,lowerArgument) # FIXME: Double check the need for " + pi" and the minus sign
        # I think that minus sign on '-upperArgument' is good. It made the angles make sense
        # Nevermind, I think it should stay positive, it made a symmetric distribution in my sources at zero AOA
        #println(Beta)
    else
        Beta = pi
    end

    return Beta

end

function getInformation(currentPanel,currentControlPoint,controlPointPanel,i,j)

    rij, rijPlusOne = calculateDistances(currentPanel,currentControlPoint)
    Beta = calculateBeta(currentPanel,currentControlPoint,i,j)
    currentPanelSine, currentPanelCosine, currentPanelLength = calculateSineAndCosine(currentPanel)
    controlPanelSine, controlPanelCosine, controlPanelLength = calculateSineAndCosine(controlPointPanel)

    # I used atan here to make sure we get the right quadrant for our angles
    theta_controlPanel = atan(controlPanelSine,controlPanelCosine)
    theta_currentPanel = atan(currentPanelSine,controlPanelCosine)

    # Trying the old way
    #theta_controlPanel = asin(controlPanelSine)
    #theta_currentPanel = asin(currentPanelSine)

    return rij, rijPlusOne, Beta, theta_controlPanel, theta_currentPanel

end

function getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j)

    rij, rijPlusOne = calculateDistances(currentPanel,currentControlPoint)
    Beta = calculateBeta(currentPanel,currentControlPoint,i,j)
    currentPanelSine, currentPanelCosine, currentPanelLength = calculateSineAndCosine(currentPanel)
    controlPanelSine, controlPanelCosine, controlPanelLength = calculateSineAndCosine(controlPointPanel)
        

    return rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine

end


########################################
# Main code block

# Make the Panels
c = 1 # chord length (keep at 1)
dphi = 0.01 # phi spacing for cosine spacing (keep small)

x = makeXarray(c,dphi) # make the cosine-spaced x-grid
e = 2 # max camber percentage
p = 4 # max camber position (multiply by 10 to get percentage)
t = 12 # max thickness percentage
alpha = 4 # Angle of attack (Degrees)
freestream = 1 # Freestream velocity (m/s)
Thickness,ybar = createNACA(e/100,p/10,t/100,x) # Get the thickness and mean camber line
panels = createPanels(x,ybar,Thickness)
numPanels = length(panels[:,1]) - 1 # Number of panels, for easy access later
println("Number of Panels: ",numPanels)
#plotAirfoil(panels,true)

# Find the Control Points
controlPoints = zeros(numPanels,2)
for i = 1:numPanels

    controlPoints[i,:] = getControlPoint(panels[i,:],panels[i+1,:])

end

# Calculating the sines and cosines of the panel relative to the x-axis
panelSines, panelCosines, panelLengths = calculateSineAndCosineOld(panels,numPanels)

# Calculating and plotting the unit normal vectors
unitNormals = calculateUnitNormals(panels,numPanels,panelSines,panelCosines)
unitTangents = calculateUnitTangents(panels,numPanels,panelSines,panelCosines)
scaleFactor = 10 # Scaling the unit normal vectors for visual aid
#plotControlPoints(panels,controlPoints)
#plotUnitNormals(controlPoints,unitNormals,numPanels,scaleFactor)
#plotUnitTangents(controlPoints,unitTangents,numPanels,scaleFactor)

#######################################################
# Create the linear system and solve it

# Create a double loop (i and j) that creates the A matrix
A = zeros(numPanels + 1, numPanels + 1) # (row,column)
b = zeros(numPanels + 1)

for i = 1:numPanels # i = row number
    #println("Working on panel: ",i)
    currentControlPoint = controlPoints[i,:] # good
    controlPointPanel = [panels[i,1:2],panels[i+1,1:2]] # good

    # Populating the majority of the 'A' matrix, excluding the Kutta Condition and vortex element
    for j = 1:numPanels # j = column number
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        
        # Getting the necessary information for the next few calculations
        rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good

        # Calculating a single point in the 'A' array
        #Aij = log(rijPlusOne/rij)*sin(theta_controlPanel - theta_currentPanel) + Beta*cos(theta_controlPanel - theta_currentPanel)
        Aij = log(rijPlusOne/rij)*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine) + Beta*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine)
        A[i,j] = Aij
    end

    # Populating the right-most column of the 'A' matrix, which deals with the vortex element
    for j = 1:numPanels # j = column number
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good

        #AiNplus1_increment = log(rijPlusOne/rij)*cos(theta_controlPanel - theta_currentPanel) - Beta*sin(theta_controlPanel - theta_currentPanel)
        AiNplus1_increment = log(rijPlusOne/rij)*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine) - Beta*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine)
        A[i,numPanels + 1] = A[i,numPanels + 1] + AiNplus1_increment
    end
    
    # Populating the bottom row of the 'A' matrix, which deals with the Kutta Condition and the sources
    # FIXME: I think there may be a problem in here
    if i == 1 || i == numPanels
        for j = 1:numPanels # j = column number

        currentPanel = [panels[j,1:2],panels[j+1,1:2]]

        rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good
        if j < numPanels/2
            controlPanelSine = -controlPanelSine
        end
        ANplus1j_increment = Beta*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine) - log(rijPlusOne/rij)*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine)
        A[numPanels + 1,j] = A[numPanels + 1,j] + ANplus1j_increment # Double check
        end

    end

    # Populating the bottom right corner of the 'A' matrix, which deals with the Kutta Condition and the vortex strength
    # FIXME: I think there may be a problem in here
    # FIXME: THERE MAY BE A PROBLEM IN HERE! For the symmetric case, shouldn't the signs be switched to cancel out? The sine functions give the same value each time through the loop
    if i == 1 || i == numPanels
        for j = 1:numPanels
            currentPanel = [panels[j,1:2],panels[j+1,1:2]]

            rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good
            if j < numPanels/2
                controlPanelSine = -controlPanelSine
            end
            ANplus1Nplus1_increment = Beta*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine) + log(rijPlusOne/rij)*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine)

            A[numPanels + 1, numPanels + 1] = A[numPanels + 1, numPanels + 1] + ANplus1Nplus1_increment
        end

    end

    # Populating the first 'N' spots in the 'b' matrix, excluding the final 'N+1' term
    # This is almost certainly correct
    for j = 1:numPanels
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good
        #println("i = ",i,", j = ",j,". controlPanelCosine = ",controlPanelCosine)
        if j < numPanels/2
            controlPanelSine = -controlPanelSine
            controlPanelCosine = -controlPanelCosine
        end

        #b[i] = 2*pi*freestream*sin(theta_controlPanel - alpha*pi/180)
        b[i] = 2*pi*freestream*(controlPanelSine*cos(alpha*pi/180) - controlPanelCosine*sin(alpha*pi/180))

    end

    # Populating the 'N+1' spot in the 'b' array
    for j = 1:numPanels
        if i == 1 || i == numPanels
            currentPanel = [panels[j,1:2],panels[j+1,1:2]]
            rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good      

            if j < numPanels/2
                controlPanelSine = -controlPanelSine
                controlPanelCosine = -controlPanelCosine
            end

            bNplus1_increment = -2*pi*freestream*(controlPanelCosine*cos(alpha*pi/180) + controlPanelSine*sin(alpha*pi/180))
            b[numPanels + 1] = b[numPanels + 1] + bNplus1_increment

        end
    end

end

# println(A) # The last corner before the N + 1 row and column should be pi, right?
#println(' ')
#println(b)

# Solving the linear system
println("Solving Linear System")
solution = A\b
println("Solution Found!")
#println(solution)
sources = solution[1:numPanels]
vortex = solution[numPanels + 1]
# vortex = 0 # If circulation is known, set it here to test for errors


#############################################
# Calculating Tangential Velocities and Pressure coefficient

# Now make a double loop to get all of the tangential velocities so we can calculate the
# pressure coefficients

tangentialVelocity = zeros(numPanels)
pressureCoefficient = zeros(numPanels)

for i = 1:numPanels

    currentControlPoint = controlPoints[i,:] # good
    controlPointPanel = [panels[i,1:2],panels[i+1,1:2]] # good    

    # Calculate the angle of the control panel
    panelSine, panelCosine, panelLength = calculateSineAndCosine(controlPointPanel)
    #println(panelSine)
    if i < numPanels/2
        panelSine = -panelSine
    end
    # I used atan here to make sure we get the right quadrant for our angles
    #theta_controlPanel = atan(panelSine,panelCosine)
    #println("Panel ",i,". Angle = ",theta_controlPanel)

    # Freestream component
    tangentialVelocity[i] = tangentialVelocity[i] + freestream*(panelCosine*cos(alpha*pi/180) + panelSine*sin(alpha*pi/180)) # good
    # println(tangentialVelocity[i])
    for j = 1:numPanels
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        rij, rijPlusOne, Beta, controlPanelSine, controlPanelCosine, currentPanelSine, currentPanelCosine = getInformation2(currentPanel,currentControlPoint,controlPointPanel,i,j) # good
        # if i < numPanels/2
        #     controlPanelSine = -controlPanelSine
        # end
        # if i < numPanels/2
        #     currentPanelSine = -currentPanelSine
        # end
        # Sources component
        tangentialVelocity[i] = tangentialVelocity[i] + 1/(2*pi)*sources[j]*(Beta*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine) - log(rijPlusOne/rij)*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine))

        # Vortex component
        tangentialVelocity[i] = tangentialVelocity[i] + vortex/(2*pi)*(Beta*(controlPanelCosine*currentPanelCosine + controlPanelSine*currentPanelSine) + log(rijPlusOne/rij)*(controlPanelSine*currentPanelCosine - controlPanelCosine*currentPanelSine))

    end

end

#println(tangentialVelocity)

pressureCoefficient = 1 .- (tangentialVelocity./freestream).^2
cp = pressureCoefficient # WARNING: Pass-by-reference problem may occur if you change pressureCoefficient later

#println(pressureCoefficient)

# Make a plottable x-axis
plotAirfoil(panels,false)
plot(controlPoints[:,1],-pressureCoefficient,color = "green",label = "Panel Method")
ylim(-1.5,2)
title("Pressure Coefficient Distribution")
xlabel("x")
ylabel("-Cp")
grid("on")
legend()

# Calculating the Lift using surface integration
cl = zeros(numPanels) # Lift Coefficients
cd = zeros(numPanels) # Drag Coefficients
gamma = zeros(numPanels) # Circulation Values
cm = zeros(numPanels,3) # Moment Coefficient
ac = [0.25,0] # Aerodynamic Center
for i = 1:numPanels

    currentPanel = [panels[i,1:2],panels[i+1,1:2]]
    #println(controlPoints[i])

    panelSine, panelCosine, panelLength = calculateSineAndCosine(currentPanel)

    # panelSine needs a correction here
    if i < numPanels/2
        panelSine = -panelSine
    end

    cl[i] = cp[i] * panelLength * panelCosine
    cd[i] = cp[i] * panelLength * panelSine
    gamma[i] = vortex*panelLength
    distance = [controlPoints[i,1] - ac[1], controlPoints[i,2] - ac[2],0]
    # println(distance)
    force = [-cd[i],-cl[i],0]
    cm[i,:] = cross(force,distance)
    # println(cm[i,:])

end

CL_integration = -sum(cl)
CD_integration = -sum(cd)
Circulation = sum(gamma)
CM_integration = sum(cm) # Only in the z-direction so this type of sum works. Investigate further if non-2D flow

# Calculating the lift using Kutta-Joukowski
# The Kutta-Joukowski Theorem is L' = density(Velocity Vector x Circulation Vector)
# L' = Lift per unit span
Velocity = [freestream*cos(alpha*pi/180),freestream*sin(alpha*pi/180),0]
Circulation = [0,0,Circulation]
density = 1 # Density of freestream
S = 1 # Unit planform area
Force = density.*cross(Velocity,Circulation)
L = density*norm(Velocity)*norm(Circulation)
dynamicPressure = .5 * density * freestream^2
forceCoefficients = Force/(dynamicPressure * S)
CL = L/(dynamicPressure * S)
#println("Force Coefficients = ",forceCoefficients)
#println("CL = ",CL*2)
CD_KuttaJouskowski = forceCoefficients[1]
CL_KuttaJouskowski = -forceCoefficients[2]

println("CL_integration: ",CL_integration)
println("CL_KuttaJouskowski: ",CL_KuttaJouskowski)
println("CD_integration: ",CD_integration)
println("CD_KuttaJouskowski: ",CD_KuttaJouskowski)
println("CM_integration: ",CM_integration)

#println(controlPoints[:,1])
println(-pressureCoefficient)
#return A # for the N x N section, each value should be symmetric about the matrix center
#return solution

# Functions that didn't make the cut because they weren't necessary
# Function (Don't use): Rotate all the coordinates into the star coordinate system for a single panel
# Function (Don't use): Calculate the induced velocity caused by the sources of one panel on another panel
# Function (Don't use): Calculate the induced velocity caused by the vortices of one panel on another panel
# Function (Don't use): Rotate velocities from the star coordinate system back into the normal coordinates