# include("Desktop/School/Aerodynamics/Homework 2/Homework 2-1 general.jl")

using Pkg
Pkg.add("PyPlot")
Pkg.add("FLOWMath")
using PyPlot
using FLOWMath

include("taf.jl")

function getCoefficients(e,p,t,x,alpha,numCoefficients)
    e = e/100
    p = p/10
    t = t/100

    alpha = alpha # Example Fourier Coefficients based on AOA = 2 Degrees

    T, ybar, dTdx, dybardx = naca4(e, p, t, x)

    # # Plot the airfoil
    # figure()
    # plot(x,ybar+0.5*T)
    # plot(x,ybar-0.5*T)
    # ylim(-.3,.3)

    integrand = zeros(1,length(phi)) # to hold the integral values before they are summed

    ###########################################
    # Solving for A0
    coefficients = zeros(numCoefficients,1)
    for i = 1:length(phi)
        integrand[i] = dybardx[i]*dphi
    end

    coefficients[1] = alpha - 1/pi * sum(integrand)

    ############################################
    # Solving for the other Fourier Coefficients
    
    for n = 2:length(coefficients)
        for i = 1:length(phi)
            integrand[i] = dybardx[i] * cos((n-1) * phi[i]) * dphi
        end

        coefficients[n] = 2/pi * sum(integrand)

    end

    println("A0 = ",coefficients[1]) # .0304
    println("A1 = ",coefficients[2]) # 0.08149
    println("A2 = ",coefficients[3]) # 0.01386
    println("A3 = ",coefficients[4]) # 0.00277
    println("A4 = ",coefficients[5]) # -0.0021

    return transpose(coefficients), T, ybar, dTdx, dybardx

end

function gammaVals(Vinf,phi,coefficients)

    gamma = zeros(length(phi))
    A0 = coefficients[1]

    gamma = 2*Vinf*A0*(1 .+ cos.(phi))./sin.(phi)

    for i = 1:length(coefficients)-1
        gamma += 2*Vinf*coefficients[i+1]*sin.(i.*phi)
    end

    return gamma

end

##############################################
# Plotting CL and Cmac

e = 2
p = 4
t = 12
alpha = 2*pi/180
numCoefficients = 5
coefficients = zeros(length(alpha),numCoefficients)
CL = zeros(length(alpha),1)
M = zeros(length(alpha),1)
c = 1
dphi = 0.01
phi = dphi:dphi:pi
x = zeros(length(phi),1)
for i = 1:length(phi)
    x[i] = c/2 * (1-cos(phi[i]))
end
rho = 1
V = 1
c = 1

coefficients, T, ybar, dTdx, dybardx = getCoefficients(e,p,t,x,alpha,numCoefficients) # checked

# Lift Coefficient
CL = 2*pi*(coefficients[1] + 0.5 * coefficients[2])

# Moment Coefficient about quarter-chord line
M = -pi/4 * (coefficients[2] - coefficients[3])

# Pressure Coefficient
qofx = V * dTdx # checked
println(dTdx)
gamma = gammaVals(V,phi,coefficients) # checked
figure()
plot(phi,gamma)
    
# Thickness velocities
ut = similar(x)
integrand = similar(x)
for i = 1:length(x)
    integrand = qofx./(x[i] .- x)
    integrand[i] = 0.0  # remove singularity
    ut[i] = 1.0/(2*pi)*trapz(x, integrand)
    #println(integrand)
end
vt_top = qofx/2 # good
vt_bottom = -qofx/2 # good

# Camber velocities
uc_top = gamma./2 # good
uc_bottom = -gamma./2 # good
vc = -V.*(alpha .- dybardx)
    
Velsqr_top = (V*cos(alpha) .+ ut .+ uc_top).^2 .+ (V*sin(alpha) .+ vt_top .+ vc).^2
Cp_top = 1 .- Velsqr_top./(V^2)

Velsqr_bottom = (V*cos(alpha) .+ ut .+ uc_bottom).^2 .+ (V*sin(alpha) .+ vt_bottom .+ vc).^2
Cp_bottom = 1 .- Velsqr_bottom./(V^2)

figure()
plot(x,ut)
title("ut")
ylabel("ut")
xlabel("x")

figure()
plot(x,uc_top)
title("uc_top")
ylabel("uc_top")
xlabel("x")

figure()
plot(x,uc_bottom)
title("uc_bottom")
ylabel("uc_bottom")
xlabel("x")

figure()
plot(x,vt_top)
title("vt_top")
ylabel("vt_top")
xlabel("x")

figure()
plot(x,vt_bottom)
title("vt_bottom")
ylabel("vt_bottom")
xlabel("x")

figure()
plot(x,vc)
title("vc")
ylabel("vc")
xlabel("x")

figure()
plot(x,qofx)
title("q")
ylabel("q")
xlabel("x")

figure()
plot(x,-Cp_top,label = "Top")
plot(x,-Cp_bottom,label = "Bottom")
title("Pressure Coefficient at AOA = 4 Degrees")
xlabel("x")
ylabel("Cp")
ylim(-1.0,2.0)
grid("on")
legend()