using Plots, LinearAlgebra

# Constants definitions
#====================================================================================#
wc = 50;

# Fucntion definitions
#====================================================================================#
# DQ Transformation
function parks_transf(Vabc, Iabc, theta_t)
    cols = size(Vabc, 2)
    Vdq0 = zeros(3, cols)
    Idq0 = zeros(3, cols)
    for i in 1:cols
        T_p = (2/3) * [cos(theta_t[i])       cos(theta_t[i] - (2 * pi / 3))   cos(theta_t[i] + (2 * pi / 3));
                       sin(theta_t[i])       sin(theta_t[i] - (2 * pi / 3))   sin(theta_t[i] + (2 * pi / 3));
                       1/2                   1/2                             1/2]
        Vdq0[:, i] = T_p * Vabc[:, i]
        Idq0[:, i] = T_p * Iabc[:, i]
    end
    return Vdq0, Idq0
end


# Inegrator step fucntion with Forward Euler
function intgstep_fwdeuler(f, y0, t0, tn, h)
    times = t0:h:tn
    y = zeros(length(times))
    y[1] = y0

    for i in 1:length(times)-1
        y[i+1] = y[i] + h * f(y[i])
    end

    return times, y
end

# Phase-locked loop with Forward Euler
function PLL_fwdeuler(theta_t, zeta0, theta_pll0, t0, tn, h)
    
    Kp_pll = 0.2
    Ki_pll = 5
    ωb = 2*pi*60
    ω0 = 1
    ωDQ = 0.98
     
    n_steps = Int((tn - t0) / h) + 1
    times = t0:h:tn
    zeta = zeros(n_steps)
    theta_pll = zeros(n_steps)
    omega_pll = zeros(n_steps)
    
        zeta[1] = zeta0
        theta_pll[1] = theta_pll0
    
        for i in 1:n_steps-1
            omega_pll[i] = Kp_pll * (theta_t - theta_pll[i]) + Ki_pll * zeta[i]
    
            # Forward Euler steps
            zeta[i+1] = zeta[i] + h * (theta_t - theta_pll[i])
            theta_pll[i+1] = theta_pll[i] + h * (omega_pll[i] + w0 - wDQ) * wb
        end
        
        return times, zeta, theta_pll, omega_pll
    end

# MAIN
#====================================================================================#
# %% Metered values at the terminal/PCC
i = 1
Vabc = zeros(3, 11) 
theta_t = zeros(11) 
for t in 0.0:0.001:0.01 
    Vabc[:, i] = [1.0 * sin(2 * pi * 60 * t);            # Phase A
                  1.0 * sin(2 * pi * 60 * t - 2*pi/3);   # Phase B
                  1.0 * sin(2 * pi * 60 * t + 2*pi/3)]   # Phase C
    theta_t[i] = 2 * pi * 60 * t + 0.25
    i += 1
end

Z = 10 
Iabc = Vabc / Z 

# %% Dynamic Phase Angle
  # Assuming a 60 Hz system

#=
    Vrms = (1/sqrt(2)) * (12.47 * 1000)
    Vabc = Vrms .* ones(3, 11)
    Z = 10
    Iabc = Vabc / Z
    theta_t = (pi / 4) .* ones(11) # (theta_t : angle between DQ frame and Vt)
=#
# Parks Transformation
Vdq0, Idq0 = parks_transf(Vabc, Iabc, theta_t)

# All these values are time variants (but the real )
Vt_d = Vdq0[1,:]
Vt_q = Vdq0[2,:]
It_d = Idq0[1,:]
It_q = Idq0[2,:]

Vt= Vt_d + im*Vt_q
It= It_d + im*It_q

# %% Power Calculation and Filtering

p = Vt_d.*It_d + Vt_q.*It_q
q = Vt_q.*It_d + Vt_d.*It_q
p_hat = zeros(length(p))
q_hat = zeros(length(q))
    for i in 1:length(p)
        f(p_h) = wc*(p[i] - p_h)
        g(q_h) = wc*(q[i] - q_h)
        times, p_h = intgstep_fwdeuler(f, 0, 0, 7, 0.01)
        times, q_h = intgstep_fwdeuler(g, 0, 0, 7, 0.01)
        p_hat[i] = p_h[end]
        q_hat[i] = q_h[end]
        #plot(times, p_h)
    end

# Phase locked loop (PLL)

theta_t = atan.(Vt_q./Vt_d)
theta_pll = zeros(length(theta_t))
omega_pll = zeros(length(theta_t))
for i in 1:length(theta_t )
    theta = theta_t[i]
    t, zeta, thetapll, omegapll = PLL_fwdeuler(0.648, 0, 0, 0, 10, 0.001)
    theta_pll[i] = thetapll[end]
    omega_pll[i] = omegapll[end]
    #plot(times, theta_pll, label="theta_pll", title="Theta_pll and Omega_pll over Time", xlabel="Time", ylabel="Value")
    #plot!(times, omega_pll, label="omega_pll")
end

# Active Power and Droop Block
P0 = 0.5
q0 = 0.1
v0= 1
mp = 100 
mq = 0.05

P_star = P0 .- mp.*omega_pll
v_star = v0 .- mq.*(q_hat .- q0)

