using Plots, LinearAlgebra

function pll_simulation(kp::Float64, ki, wDQ, theta_t)
    #kp = 0.8
    ki = 0.5
    wDQ = 1
    theta_t = 10*(pi/180) 

    w0 = 1.0                   # pu
    wb = 2*pi* 60.0            # base ang freq

    A = wb*[-kp ki; -1/wb 0]
    B = wb*[kp 1 -1; 1/wb 0 0]
    C = [1 0; -kp ki]
    D = [0 0 0; kp 0 0]

    U = [theta_t w0 wDQ]'

    x_dot = a*x + b*U
    y = ..

    return x_dot, y

end

    dt = 0.0001
    T = 5.0



    # AX + BU = 0 for initial conditions of X

    x0 = A \ (B * U)
    y0 = C*x0 + D*U

    xF = zeros(length(x0), Int(T / dt) + 1)
    tF = zeros(length(x0), Int(T / dt) + 1)
    yF = zeros(length(y0), Int(T / dt) +1 )

    # Set the initial condition in the first column
    xF[:, 1] .= x0
    tF[1] = 0.0
    yF[:, 1] .= y0

    for k in 1:Int(round(T / dt))
        tF[:, k + 1] .= tF[:, k] .+ dt
        xF[:, k + 1] .= xF[:, k] + dt * (A * xF[:, k] + B * U)
        yF[:, k + 1] .= yF[:, k] + dt * (C * yF[:, k] + D * U)
    end

    #return tF[2, :], xF[1, :], xF[2, :]
#end

plot(tF[2, :], yF[1, :], label="theta_pll")
plot!(tF[2, :], yF[2, :], label="omega_pll")

# the function pll_simulation
#tf, theta_pll, omega_pll = pll_simulation(0.8, 0.5, 1, π / 4)

#plot(t, theta_pll, label="theta_pll")
#plot!(t, omega_pll, label="omega_pll")