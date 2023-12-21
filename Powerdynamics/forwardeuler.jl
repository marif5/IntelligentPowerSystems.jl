wc = 50
h = 0.01
y0 = 0
t0 = 0
tn = 7

f(p) = wc*(200 - p)
times = t0:h:tn
y = zeros(length(times))
y[1] = y0

    for i in 1:length(times)-1
        y[i+1] = y[i] + h * f(y[i])
        
    end
plot(times, y)


#==================================================#

function solve_pll(theta_t, zeta0, theta_pll0, t0, tn, h)
    
Kp_pll = 0.2
Ki_pll = 5
ωb = 2*pi*60
ω0 = 1
ωDQ = 0.95
 
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

t0 = 0 # Start time
tn = 10 # End time
h = 0.01 # Time step
theta_t = 0.76 

times, zeta, theta_pll, omega_pll = solve_pll(theta_t, 0, 0, t0, tn, h)
plot(times, theta_pll, label="theta_pll", title="Theta_pll and Omega_pll over Time", xlabel="Time", ylabel="Value")
plot!(times, omega_pll, label="omega_pll")

#==================================================#

function VI_conttollers(theta_t, zeta0, theta_pll0, t0, tn, h)
    
    Kp_v = 0.2
    Ki_v = 5
    Kf_v = 
    Cf = 
    ωb = 2*pi*60
    ω0 = 1
    ωDQ = 0.95
     
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
    
    t0 = 0 # Start time
    tn = 10 # End time
    h = 0.01 # Time step
    theta_t = 0.76 
    
    times, zeta, theta_pll, omega_pll = solve_pll(theta_t, 0, 0, t0, tn, h)
    plot(times, theta_pll, label="theta_pll", title="Theta_pll and Omega_pll over Time", xlabel="Time", ylabel="Value")
    plot!(times, omega_pll, label="omega_pll")

    #=============================#

fucntion inverter_der(X, xline, omega_c, Kp_pll, Ki_pll, Kp_i, omega_0, omega_DQ, omega_b, p0, q0, v0, mp, mq)
    p_hat       = X(1)
    q_hat       = X(2)
    zeta        = X(3)
    theta_pll   = X(4)
    delta       = X(5)
    phi_d       = X(6)
    gamma       = X(7)
    vs          = X(8)
    Id_s        = X(9)
    Iq_s        = X(10)
    Vd_t        = X(11)
    Vq_t        = X(12)

    iq_t = (1-Vd_t)/xline
    id_t = Vq_t/xline

    p = Vd_t*iq_t + Vq_t*iq_t
    q = Vq_t*iq_t + Vd_t*iq_t

    p_hat_dot = omega_c*(p-p_hat)
    q_hat_dot = omega_c*(q-q_hat)

    theta_t = atan(Vq_t./Vd_t)
    
    zeta_dot = theta_t - theta_pll
    omega_pll = Kp_pll*(theta_t-theta_pll)+Ki_pll*zeta
    theta_pll_dot = (omega_pll + omega_0 - omega_DQ)*omega_b

    p_star = p0 - mp*omega_pll
    q_star = q0 - mq*(q_hat - q0)

end