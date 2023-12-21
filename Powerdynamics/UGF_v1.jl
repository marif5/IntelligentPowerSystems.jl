using Plots, LinearAlgebra, JuMP, Ipopt

# Constants definitions
#====================================================================================#
wc = 50;

# Fucntion definitions
#====================================================================================#
# DQ Transformation
function parks_transf(Vabc, theta_t)
    cols = size(Vabc, 2)
    Vdq0 = zeros(3, cols)
    for i in 1:cols
        T_p = (2/3) * [cos(theta_t[i])       cos(theta_t[i] - (2 * pi / 3))   cos(theta_t[i] + (2 * pi / 3));
                       sin(theta_t[i])       sin(theta_t[i] - (2 * pi / 3))   sin(theta_t[i] + (2 * pi / 3));
                       1/2                   1/2                             1/2]
        Vdq0[:, i] = T_p * Vabc[:, i]
    end
    return Vdq0
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

Vdq0 = parks_transf(Vabc, theta_t)
Vdt = Vdq0[1,1]
Vqt = Vdq0[2,1]

#=========== Dr.Sam's Approach to Fowrward Euler =============#

function inverter_der(X, xline, V_inf0, PLL_Coefficient, P_filter_droop_Coefficient, IV_Controller_Coefficient)
    #stating coefficents:
        Kp_pll      = PLL_Coefficient[1]
        Ki_pll      = PLL_Coefficient[2]
        omega_0     = PLL_Coefficient[3]
        omega_DQ    = PLL_Coefficient[4]
        omega_b     = PLL_Coefficient[5]

        omega_c     = P_filter_droop_Coefficient[1]
        Kp_i        = P_filter_droop_Coefficient[2]
        p0          = P_filter_droop_Coefficient[3]
        q0          = P_filter_droop_Coefficient[4]
        v0          = P_filter_droop_Coefficient[5]
        mp          = P_filter_droop_Coefficient[6]
        mq          = P_filter_droop_Coefficient[7]
    
        Kp_vc       = IV_Controller_Coefficient[1]
        Ki_vc       = IV_Controller_Coefficient[2]
        Kf_vc       = IV_Controller_Coefficient[3]
        Kp_ic       = IV_Controller_Coefficient[4]
        Ki_ic       = IV_Controller_Coefficient[5]
        Kf_ic       = IV_Controller_Coefficient[6]
        Cf          = IV_Controller_Coefficient[7]
        Lf          = IV_Controller_Coefficient[8]
        Kv_i        = IV_Controller_Coefficient[9]

# Forward Euler formation    
    h = 0.02    
    times = 0.02:h:0.04
    for i in 1:length(times)
        #stating states:
        p_hat       = X[1]      #powerfilter
        q_hat       = X[2]      #powerfilter
        zeta        = X[3]      #PLL
        theta_pll   = X[4]      #PLL
        delta       = X[5]      #Active Power Control
        phi_d       = X[6]      #IV Contollers
        gamma_d     = X[7]      #IV Contollers
        Vs_abs      = X[8]      #IV Contollers
        Id_s        = X[9]      #IV Contollers - ??
        Iq_s        = X[10]     #IV Contollers - ??
        Vd_t        = X[11]     #IV Contollers - Non Zero
        Vq_t        = X[12]     #IV Contollers - Non Zero

      
        #state equations
        Iq_t = (V_inf0-Vd_t)/xline                                  # Non-state variable - valid   
        Id_t = Vq_t/xline                                           # Non-state variable - valid

        p           = Vd_t*Id_t + Vq_t*Iq_t                         # Non-state variable - valid
        q           = Vq_t*Id_t - Vd_t*Iq_t                         # Non-state variable - valid
        p_hat_dot   = omega_c*(p-p_hat)                             # substate variable - valid
        q_hat_dot   = omega_c*(q-q_hat)                             # substate variable - valid

        theta_t     = atan(Vq_t/Vd_t)                               # Non-state variable - valid
        
        zeta_dot       = theta_t - theta_pll                        # Non-state variable - valid
        omega_pll      = Kp_pll*(theta_t-theta_pll) + Ki_pll*zeta   # Non-state variable - valid
        theta_pll_dot  = (omega_pll + omega_0 - omega_DQ)*omega_b   # substate variable - valid

        p_star      = p0 - mp*omega_pll                             # Non-state variable - valid
        delta_dot   = Kp_i*(p_star - p_hat)                         # Non-state variable - valid
        v_star      = v0 - mq*(q_hat - q0)                          # Non-state variable - valid

        phi_d_dot   = v_star - Vd_t                                 # substate variable - valid
        Id_s_hat    = Kp_vc*phi_d_dot - Ki_vc*phi_d  - Kf_vc*Id_t - (omega_pll - omega_0)*Cf*Vq_t           # Non-state variable - valid

        gamma_d_dot = Id_s_hat - Id_s                               # substate variable - valid
        Vd_s_hat    = Kp_ic*gamma_d_dot - Ki_ic*gamma_d  - Kf_ic*Vd_t - (omega_pll - omega_0)*Lf*Iq_s       # Non-state variable - valid

        Vt_mag      = sqrt((Vd_t)^2 + (Vq_t)^2)
        Vs_abs_dot  = Kv_i*(v_star - Vt_mag)                        # Non-state variable - valid
        
        Vd_s    =   Vs_abs*cos(delta)
        Vq_s    =   Vs_abs*sin(delta)
        #Vd_s    =   0
        #Vq_s    =   0
        Id_s_dot    = (omega_b*(Vd_s - Vd_t) + (omega_pll + omega_0)*omega_b*Iq_s)/Lf;      # Non-state variable - valid
        Iq_s_dot    = (omega_b*(Vq_s - Vq_t) - (omega_pll + omega_0)*omega_b*Id_s)/Lf;      # Non-state variable - valid
        Vd_t_dot    = (omega_b*(Id_s - Id_t) + (omega_pll + omega_0)*omega_b*Vq_t)/Cf;      # Non-state variable - valid
        Vq_t_dot    = (omega_b*(Iq_s - Iq_t) - (omega_pll + omega_0)*omega_b*Vd_t)/Cf;      # Non-state variable - valid
        

        fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
        X  = X + h* fX

    end    

    return times, X       
end

function init_inv(p0, q0, v0, omega_0, Lf, Cf, xline)
    model = Model(Ipopt.Optimizer)
    @variable(model, (-pi/4) <= theta_t_0 <= (pi/4))
    @variable(model, 0.5 <= v_inf <= 1.5)
    @variable(model, (-pi/4) <= delta <= (pi/4))
    @variable(model, 0.5 <= Vs0 <= 1.5)

    # constraints
    @constraint(model, p0 == (v0*v_inf/xline)*sin(theta_t_0))
    @constraint(model, q0 == ((v0*v_inf*cos(theta_t_0) - (v_inf^2))/xline))
    @constraint(model, p0 == (v0*Vs0/(omega_0*Lf))*sin(delta-theta_t_0))
    @constraint(model, q0 == (v0*Vs0*cos(delta-theta_t_0) - (Vs0^2))/(omega_0*Lf) +(v0^2)*(Cf*omega_0))

    @objective(model, Min, 1)
    optimize!(model)

    # make sure the solution is valid
    println(termination_status(model))
    println(Int(termination_status(model)))
    if (Int(termination_status(model)) != 4)
        @assert 1 == 2
    end
    return value(v_inf),value(theta_t_0), value(Vs0), value(delta)
end

#=== Main ==#



line_X       = 0.2;
omega_c     = 50;
Kp_pll      = 0.2;
Ki_pll      = 5.0;
Kp_i        = 0.3;
omega_0     = 1;
omega_DQ    = 0.98;
omega_b     = 2*pi*60;
p0          = 0.5;
q0          = 0.1;
v0          = 1;
mp          = 100;
mq          = 0.05;
Kp_vc       = 1;
Ki_vc       = 2 ;
Kf_vc       = 1;
Kp_ic       = 1;
Ki_ic       = 2 ;
Kf_ic       = 0;
Cf          = 0.074;
Lf          = 0.08;
Kv_i        = Kp_i; #not given in the paper

V_inf0, theta_t0, Vs0, delta0 = init_inv(p0, q0, v0, omega_0, Lf, Cf, line_X)
# %%
    X = zeros(12);
    X[1] = p0
    X[2] = q0
    X[3] = 0                    #zeta
    X[4] = theta_t0             #theta_pll
    X[5] = delta0               #delta
    X[6] = 0                    #phi_d
    X[7] = 0                    #gamma_d
    X[8] = 1                    #Vs_abs
    X[11] = 1*cos(theta_t0)     #Vd_t
    X[12] = 1*sin(theta_t0)     #Vq_t

    Iq_t = (V_inf0-X[11])/line_X ;     
    Id_t = X[12]/line_X ;              
    omega_pll_0 = Kp_pll*X[3];
    X[9] = -1*(omega_pll_0 + omega_0)*X[12]*Cf + Id_t    #Id_s
    X[10] = ((omega_pll_0 + omega_0)*X[11]*Cf) + Iq_t    #Iq_s

PLL_coeff = [Kp_pll, Ki_pll, omega_0, omega_DQ, omega_b]
P_droop_Coeff = [omega_c, Kp_i, p0, q0, v0, mp, mq]
IVControl_coeff = [Kp_vc, Ki_vc, Kf_vc, Kp_ic, Ki_ic, Kf_ic, Cf, Lf, Kv_i]

# Calling the function
times, states = inverter_der(X, line_X, V_inf0, PLL_coeff, P_droop_Coeff, IVControl_coeff)
states