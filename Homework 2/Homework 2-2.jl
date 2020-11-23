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
function plotAirfoil(panels)

    # Plot the airfoil
    figure()
    scatter(panels[:,1],panels[:,2],label = "Panel Tips")
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
    legend()

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
    deltaY = currentPanel[2][2,1] - currentPanel[1][2,1]
    length = sqrt(deltaX^2 + deltaY^2)
    panelLength = length
    panelSine = deltaY/length
    panelCosine = deltaX/length

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

    # plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    deltaX = currentPanel[2][1,1] - currentControlPoint[1]
    deltaY = currentPanel[2][2,1] - currentControlPoint[2]
    rijPlusOne = sqrt(deltaX^2 + deltaY^2)

    # plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    return rij, rijPlusOne

end

# Function (Use): Make sure that the distances calculated are correct
function plotDistances(currentPanel,currentControlPoint,deltaX,deltaY)

    X = [currentControlPoint[1],currentControlPoint[1] + deltaX]
    Y = [currentControlPoint[2],currentControlPoint[2] + deltaY]

    plot(X,Y,color = "black")

end

# Function (Might Use): Calculate "r" vectors
# function calculateRVectors(panels,numPanels,controlPoints)

#     for i = 1:numPanels
#         currentControlPoint = controlPoints[i,:]
#         for j = 1:numPanels
#             if i != j
#                 currentPanel = [panels[j,1:2],panels[j+1,1:2]]
#                 rij, rijPlusOne = calculateDistances(currentPanel,currentControlPoint)
#             end
#         end
#     end

#     return rij, rijPlusOne

# end

# Function (Use): Compute the angle Beta
function calculateBeta(currentPanel,currentControlPoint,i,j)

    # currentPanel = [X1 Y1; X2 Y2]
    # currentControlPoint = [Xc Yc]

    xj = currentPanel[1][1,1]
    xjPlus1 = currentPanel[2][1,1]

    yj = currentPanel[1][2,1]
    yjPlus1 = currentPanel[2][2,1]

    xi = currentControlPoint[1]
    yi = currentControlPoint[2]

    if i != j
        upperArgument = (xj - xi) * (yjPlus1 - yi) - (yj - yi) * (xjPlus1 - xi)
        lowerArgument = (xj - xi) * (xjPlus1 - xi) - (yj - yi) * (yjPlus1 - yi)
        Beta = atan(upperArgument,lowerArgument)
    else
        Beta = pi
    end

    return Beta

end


########################################
# Main code block

# Make the Panels
c = 1 # chord length (keep at 1)
dphi = 0.05 # phi spacing for cosine spacing (keep small)
x = makeXarray(c,dphi) # make the cosine-spaced x-grid
e = 2 # max camber percentage
p = 4 # max camber position (multiply by 10 to get percentage)
t = 12 # max thickness percentage
alpha = 2 # Angle of attack (Degrees)
freestream = 1 # Freestream velocity (m/s)
Thickness,ybar = createNACA(e/100,p/10,t/100,x) # Get the thickness and mean camber line
panels = createPanels(x,ybar,Thickness)
numPanels = length(panels[:,1]) - 1 # Number of panels, for easy access later
#plotAirfoil(panels)

# Find the Control Points
controlPoints = zeros(numPanels,2)
for i = 1:numPanels

    controlPoints[i,:] = getControlPoint(panels[i,:],panels[i+1,:])

end
plotControlPoints(panels,controlPoints)

# Calculating the sines and cosines of the panel relative to the x-axis
panelSines, panelCosines, panelLengths = calculateSineAndCosine(panels,numPanels)

# Calculating and plotting the unit normal vectors
unitNormals = calculateUnitNormals(panels,numPanels,panelSines,panelCosines)
unitTangents = calculateUnitTangents(panels,numPanels,panelSines,panelCosines)
scaleFactor = 10 # Scaling the unit normal vectors for visual aid
plotUnitNormals(controlPoints,unitNormals,numPanels,scaleFactor)
plotUnitTangents(controlPoints,unitTangents,numPanels,scaleFactor)

# Create a double loop (i and j) that creates the A matrix
A = zeros(numPanels + 1, numPanels + 1) # (row,column)
b = zeros(numPanels + 1)
ANplus1j = zeros(numPanels)
ANplus1Nplus1 = zeros(numPanels)
bNplus1 = zeros(numPanels)
for i = 1:numPanels
    # Select the control point that we will be using for this round
    currentControlPoint = controlPoints[i,:]
    controlPointPanel = [panels[i,1:2],panels[i+1,1:2]]
    for j = 1:numPanels
        # Select the panel that we will be using for this particular iteration
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        rij, rijPlusOne = calculateDistances(currentPanel,currentControlPoint)
        Beta = calculateBeta(currentPanel,currentControlPoint,i,j)
        currentPanelSine, currentPanelCosine, currentPanelLength = calculateSineAndCosine(currentPanel)
        controlPanelSine, controlPanelCosine, controlPanelLength = calculateSineAndCosine(controlPointPanel)

        # If this doesn't work, use eqn 2.196-8
        theta_i = asin(controlPanelSine)
        theta_j = asin(currentPanelSine)

        # Add to the 'A' matrix
        Aij = log(rijPlusOne/rij) * sin(theta_i - theta_j) + Beta * cos(theta_i - theta_j)
        AiNPlus1 = log(rijPlusOne/rij) * cos(theta_i - theta_j) - Beta * sin(theta_i - theta_j)
        A[i,j] = Aij
        # I've made edits to the next line
        A[i,numPanels + 1] = A[i,numPanels + 1] + AiNPlus1

        # Add to the 'b' matrix
        b[i] = 2*pi*freestream * sin(theta_i - alpha*pi/180)

        # Implementing the Kutta Condition at the first and last panels
        if i == 1 || i == numPanels
            # This will be a long array with information only in the first and last elements
            # The array will be summed at the end and all the zeros won't matter
            ANplus1j[i] = Beta * sin(theta_i - theta_j) - log(rijPlusOne/rij) * cos(theta_i - theta_j)
            # Same for this array
            # BUT... we will need to make an additional sum over all of the j panels
            ANplus1Nplus1[i] = Beta * cos(theta_i - theta_j) + log(rijPlusOne/rij) * sin(theta_i - theta_j)

            # Same for this array
            bNplus1[i] = -2*pi*(freestream) * cos(theta_i - alpha*pi/180)

        end

        # Adding the Kutta conidtion to the bottom row of the 'A' matrix
        A[numPanels + 1,j] = sum(ANplus1j)
        A[numPanels + 1,numPanels + 1] = sum(ANplus1Nplus1)

        # Adding the Kutta condition to the last element of the 'b' array
        b[numPanels + 1] = sum(bNplus1)

    end

end

# Solving the linear system
solution = A\b

sources = solution[1:numPanels]
vortex = solution[numPanels + 1]

# Now make a double look to get all of the tangential velocities so we can calculate the
# pressure coefficients

tangentialVelocity = zeros(numPanels)
pressureCoefficient = zeros(numPanels)
for i = 1:numPanels
    # Select the control point that we will be using for this round
    currentControlPoint = controlPoints[i,:]
    controlPointPanel = [panels[i,1:2],panels[i+1,1:2]]

    controlPanelSine, controlPanelCosine, controlPanelLength = calculateSineAndCosine(controlPointPanel)
    theta_i = asin(controlPanelSine)
    tangentialVelocity[i] = tangentialVelocity[i] + freestream * cos(theta_i - alpha*pi/180)

    #println("Tangential Velocity at panel ",i,": ",tangentialVelocity[i])

    for j = 1:numPanels

        # Select the panel that we will be using for this particular iteration
        currentPanel = [panels[j,1:2],panels[j+1,1:2]]
        rij, rijPlusOne = calculateDistances(currentPanel,currentControlPoint)
        Beta = calculateBeta(currentPanel,currentControlPoint,i,j)
        currentPanelSine, currentPanelCosine, currentPanelLength = calculateSineAndCosine(currentPanel)
        controlPanelSine, controlPanelCosine, controlPanelLength = calculateSineAndCosine(controlPointPanel)

        # If this doesn't work, use eqn 2.196-8
        theta_i = asin(controlPanelSine) # checked
        theta_j = asin(currentPanelSine) # checked


        # println("Source at panel ",j,": ",sources[j])
        tangentialVelocity[i] = tangentialVelocity[i] + 1/(2*pi) * sources[j] * ((Beta*sin(theta_i - theta_j) - log(rijPlusOne/rij)*cos(theta_i - theta_j))) # checked
        tangentialVelocity[i] = tangentialVelocity[i] + vortex/(2*pi) * (Beta*cos(theta_i - theta_j) + log(rijPlusOne/rij)*sin(theta_i - theta_j)) # checked

        #println("Tangential Velocity at panel ",i,": ",tangentialVelocity[i])

    end

    println("Tangential Velocity at panel ",i,": ",tangentialVelocity[i])
    pressureCoefficient[i] = 1 - (tangentialVelocity[i]/freestream)^2

end

#println(tangentialVelocity)


# Make a plottable x-axis
figure()
plot(controlPoints[:,1],-pressureCoefficient)

# Functions that didn't make the cut because they weren't necessary
# Function (Don't use): Rotate all the coordinates into the star coordinate system for a single panel
# Function (Don't use): Calculate the induced velocity caused by the sources of one panel on another panel
# Function (Don't use): Calculate the induced velocity caused by the vortices of one panel on another panel
# Function (Don't use): Rotate velocities from the star coordinate system back into the normal coordinates