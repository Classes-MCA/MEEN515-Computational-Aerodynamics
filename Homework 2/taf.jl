function naca4(e, p, t, x)

    T = similar(x)
    dTdx = similar(x)
    ybar = similar(x)
    dybardx = similar(x)

    @. T = 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3537*x^2 + 0.2843*x^3 - 0.1015*x^4)
    @. dTdx = 10*t*(0.2969*0.5/sqrt(x) - 0.126 - 0.3537*2*x + 0.2843*3*x^2 - 0.1015*4*x^3)

    for i = 1:length(x)
        if x[i] <= p
            ybar[i] = e/p^2 * (2*p*x[i] - x[i]^2)
            dybardx[i] = e/p^2 * (2*p - 2*x[i])
        else
            ybar[i] = e/(1-p)^2 * (1 - 2*p + 2*p*x[i] - x[i]^2)
            dybardx[i] = e/(1-p)^2 * (2*p - 2*x[i])
        end
    end

    return T, ybar, dTdx, dybardx
end