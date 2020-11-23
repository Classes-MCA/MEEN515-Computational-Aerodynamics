# include("Desktop/School/Aerodynamics/Homework 2/Homework 2-1 general.jl")

using Pkg
Pkg.add("PyPlot")
using PyPlot

include("taf.jl")

function getCoefficients(e,p,t,x,alpha,numCoefficients)
    e = e/100
    p = p/10
    t = t/100

    alpha = alpha*pi/180 # Example Fourier Coefficients based on AOA = 2 Degrees

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

function gammaVals(V,phi,coefficients)

    gamma = zeros(length(phi))
    A0 = coefficients[1]

    # Performing the summation
    for i = 1:length(gamma)

        gamma[i] = A0*(1+cos(phi[i]))/sin(phi[i])

        for n = 2:length(coefficients[1,:])

            gamma[i] = gamma[i] + coefficients[n]*sin(n*phi[i])

        end

    end

    gamma = 2*V.*gamma

    return gamma

end

##############################################
# Plotting CL and Cmac

e = 2
p = 4
t = 12
alpha = 4:4
numCoefficients = 10
coefficients = zeros(length(alpha),numCoefficients)
CL = zeros(length(alpha),1)
M = zeros(length(alpha),1)
c = 1
dphi = 0.01
phi = 0:dphi:pi
x = zeros(length(phi),1)
for i = 1:length(phi)
    x[i] = c/2 * (1-cos(phi[i]))
end
Cp_top = zeros(length(alpha),length(x))
Cp_bottom = zeros(length(alpha),length(x))
ut = zeros(length(alpha),length(x))
vt_top = zeros(length(alpha),length(x))
vt_bottom = zeros(length(alpha),length(x))

uc_top = zeros(length(alpha),length(x))
uc_bottom = zeros(length(alpha),length(x))
vc = zeros(length(alpha),length(x))
rho = 1
V = 1
c = 1

for i = 1:length(alpha)
    println("Working on angle: ",alpha[i])
    
    coefficients[i,:], T, ybar, dTdx, dybardx = getCoefficients(e,p,t,x,alpha[i],numCoefficients) # checked

    # Lift Coefficient
    CL[i] = 2*pi*(coefficients[i,1] + 0.5 * coefficients[i,2])

    # Moment Coefficient about quarter-chord line
    M[i] = -pi/4 * (coefficients[i,2] - coefficients[i,3])

    # Pressure Coefficient
    qofx = V * dTdx # checked
    gammaofphi = gammaVals(V,phi,coefficients[i,:]) # checked
    
    # Thickness velocities
    for s = 2:length(x) # s-value

        for w = 1:length(x) # x-value
            if s != w
            ut[i,s] = ut[i,s] + 1/(2*pi) * qofx[s]/(w - s) * (x[s]-x[s-1]) # See integral Dr. Ning sent out (matlab)
            end
        end

    end
    vt_top[i,:] = qofx/2 # good
    vt_bottom[i,:] = -qofx/2 # good

    # Camber velocities
    uc_top[i,:] = gammaofphi./2 # good
    uc_bottom[i,:] = -gammaofphi./2 # good
    for s = 2:length(x) # s-value

        for w = 1:length(x) # x-value
            if s != w
            vc[i,s] = vc[i,s] + 1/(2*pi) * qofx[s]/(w - s) * (x[s]-x[s-1])
            # No integral needed
            # Vc occured in the BC, when we forced it to be true
            # Vc = - Vinf *(alpha - dybardx)(eqn 2.73)
            end
        end

    end

    for j = 1:length(dybardx)
        vc[i] = -V.*(alpha[i] .- dybardx[j])
    end

    # Vmod_top    = V*(1/sqrt(1 .+ (dybardx .+ .5 * dTdx).^2))
    # Vmod_bottom = V*(1/sqrt(1 .+ (dybardx .- .5 * dTdx).^2))

    # Cp_top = 1 .- (Vmod_top./V).^2
    # Cp_bottom = 1 .- (Vmod_bottom./V).^2

    # println(Cp_bottom)

    # println("ut: ",ut[i,:])
    # println("uc_top: ", uc_top[i,:])
    # println("vt_top: ",vt_top[i,:])
    # println("vc: ",vc[i,:])
    # println("uc_bottom: ",uc_bottom[i,:])
    # println("vt_bottom: ",vt_bottom[i,:])
    
    Velsqr_top = (V*cos(alpha[i]*pi/180) .+ ut[i,:] .+ uc_top[i,:]).^2 + (V*sin(alpha[i]*pi/180) .+ vt_top[i,:] .+ vc[i,:]).^2
    Cp_top[i,:] = 1 .- Velsqr_top./(V^2)

    Velsqr_bottom = (V*cos(alpha[i]*pi/180) .+ ut[i,:] .+ uc_bottom[i,:]).^2 + (V*sin(alpha[i]*pi/180) .+ vt_bottom[i,:] .+ vc[i,:]).^2
    Cp_bottom[i,:] = 1 .- Velsqr_bottom./(V^2)

    # Vmod_top = zeros(length(dybardx))
    # Vmod_bottom = zeros(length(dybardx))
    # for j = 1:length(dybardx)
    #     Vmod_top[j]    = sqrt(Velsqr_top[j])*(1/sqrt(1 + (dybardx[j] - .5 * dTdx[j])^2))
    #     Vmod_bottom[j] = sqrt(Velsqr_bottom[j])*(1/sqrt(1 + (dybardx[j] + .5 * dTdx[j])^2))

    #     Cp_top[j] = 1 - (Vmod_top[j]./V).^2
    #     Cp_bottom[j] = 1 - (Vmod_bottom[j]./V).^2
    # end

end

figure()
scatter(alpha,CL)
plot(alpha,2*pi*(alpha .- alpha[1])*pi/180)
title("Lift Coefficient")
xlabel("Angle of Attack (Degrees)")
ylabel("Lift Coefficient")
grid("on")

#println(CL)

figure()
scatter(alpha,M)
title("Moment Coefficient about Quarter-Chord Line")
xlabel("Angle of Attack (Degrees)")
ylabel("Moment Coefficient")
grid("on")

#println(M)

figure()
plot(x,-transpose(Cp_top),label = "Top")
plot(x,-transpose(Cp_bottom),label = "Bottom")
title("Pressure Coefficient at AOA = 4 Degrees")
xlabel("x")
ylabel("Cp")
ylim(-1.0,2.0)
grid("on")
legend()