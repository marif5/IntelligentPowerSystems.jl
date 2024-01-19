# Fucntion definitions
#====================================================================================#
# DQ Transformation function

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

#====================================================================================#
# IPOT initialization of state vector function

function init_inv(p0, q0, v_t0, theta_inf, omega_0, Lf, Cf, xline)
    model = Model(Ipopt.Optimizer)
    @variable(model, (-pi/4) <= theta_t_0 <= (pi/4))
    @variable(model, (-pi/4) <= delta0 <= (pi/4))
    @variable(model, 0.5 <= v_inf <= 1.5)
    @variable(model, 0.5 <= Vs0 <= 1.5)

    # constraints
    @constraint(model, p0 == (v_t0*v_inf/xline)*sin(theta_t_0-theta_inf) )
    @constraint(model, q0 == ((-v_t0*v_inf*cos(theta_inf-theta_t_0) + (v_t0^2))/xline) )

    @constraint(model, p0 == (v_t0*Vs0/(omega_0*Lf))*sin(delta0-theta_t_0) )
    @constraint(model, q0 == (v_t0*Vs0*cos(delta0-theta_t_0) - (v_t0^2))/(omega_0*Lf) + (v_t0^2)*(Cf*omega_0) )

    @objective(model, Min, 1)
    optimize!(model)

    # make sure the solution is valid
    println(termination_status(model))
    println(Int(termination_status(model)))
    if (Int(termination_status(model)) != 4)
        @assert 1 == 2
    end
    return value(v_inf),value(theta_t_0), value(Vs0), value(delta0)
end

#====================================================================================#
# Forward Euler inverter System Model function

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
    times = 0.02:h:0.4
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
        Id_t = (Vd_t-V_inf0)/xline                                  # Non-state variable - valid   
        Iq_t = Vq_t/xline                                           # Non-state variable - valid

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
        Id_s_dot    = (omega_b*(Vd_s - Vd_t))/Lf + (omega_pll + omega_0)*omega_b*Iq_s               # Non-state variable - valid
        Iq_s_dot    = (omega_b*(Vq_s - Vq_t))/Lf - (omega_pll + omega_0)*omega_b*Id_s               # Non-state variable - valid
        Vd_t_dot    = (omega_b*(Id_s - Id_t))/Cf + (omega_pll + omega_0)*omega_b*Vq_t               # Non-state variable - valid
        Vq_t_dot    = (omega_b*(Iq_s - Iq_t))/Cf - (omega_pll + omega_0)*omega_b*Vd_t               # Non-state variable - valid
        

        fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
        X  = X + h* fX

        X1 = X2 = X3 = X4 = X5 = X6 = X6 = X7 = X8 = X9 = X10 = X11 = X12 = zeros(length(times))
        
        X1[i] = X[1];
        X2[i] = X[2];
        X3[i] = X[3];
        X4[i] = X[4];
        X5[i] = X[5];
        X6[i] = X[6];
        X7[i] = X[7];
        X8[i] = X[8];
        X9[i] = X[9];
        X10[i] = X[10];
        X11[i] = X[11];
        X12[i] = X[12];

    end    

    fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
    X  = X + h* fX
    X1 = X2 = X3 = X4 = X5 = X6 = X6 = X7 = X8 = X9 = X10 = X11 = X12 = zeros(length(times))
        X1[i] = X[1];
        X2[i] = X[2];
        X3[i] = X[3];
        X4[i] = X[4];
        X5[i] = X[5];
        X6[i] = X[6];
        X7[i] = X[7];
        X8[i] = X[8];
        X9[i] = X[9];
        X10[i] = X[10];
        X11[i] = X[11];
        X12[i] = X[12];

    Out_X = DataFrame(P_hat = X1, Q_hat = X2, zeta = X3, theta_t = X4, delta = X5, phi_d = X6, gamma_d = X7, Vs_abs = X8, Ids = X9, Iqs = X10, Vdt = X11, Vqt = X12)
    return times, X, Out_X       
end

#====================== Archived Fucntions ===================================#

function old_init_inv(p0, q0, v0, omega_0, Lf, Cf, xline)
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
#V_inf0, theta_t0, Vs0, delta0 = init_inv(p0, q0, v0, omega_0, Lf, Cf, line_X)

function inverter_dynamics(X, xline, V_inf0, PLL_Coefficient, P_filter_droop_Coefficient, IV_Controller_Coefficient)
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

    Iq_t = (V_inf0-Vd_t)/xline                                  # Non-state variable - valid   
    Id_t = Vq_t/xline                                           # Non-state variable - valid


                ##state equations
                #IQ_t = (V_inf0-Vd_t)/xline                                  # Non-state variable - valid   
                #ID_t = Vq_t/xline                                           # Non-state variable - valid
#
                #R    = [cos(theta_t0) -sin(theta_t0); sin(theta_t0) cos(theta_t0)]
                #Idq  = R\[ID_t; IQ_t]
                #Id_t = Idq[1]
                #Iq_t = Idq[2]

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
    Id_s_hat    = Kp_vc*(v_star - Vd_t) + Ki_vc*phi_d  + Kf_vc*Id_t - (omega_pll + omega_0)*Cf*Vq_t           # Non-state variable - valid

    gamma_d_dot = Id_s_hat - Id_s                               # substate variable - valid
    Vd_s_hat    = Kp_ic*(Id_s_hat - Id_s  ) + Ki_ic*gamma_d  + Kf_ic*Vd_t - (omega_pll + omega_0)*Lf*Iq_s       # Non-state variable - valid

    Vt_mag      = sqrt((Vd_t)^2 + (Vq_t)^2)
    Vs_abs_dot  = Kv_i*(v_star - Vt_mag)                        # Non-state variable - valid
    
    Vd_s    =   Vs_abs*cos(delta)
    Vq_s    =   Vs_abs*sin(delta)

    #Vd_s    =   0
    #Vq_s    =   0
    Id_s_dot    = (omega_b*(Vd_s - Vd_t))/Lf + (omega_pll + omega_0)*omega_b*Iq_s               # Non-state variable - valid
    Iq_s_dot    = (omega_b*(Vq_s - Vq_t))/Lf - (omega_pll + omega_0)*omega_b*Id_s               # Non-state variable - valid
    Vd_t_dot    = (omega_b*(Id_s - Id_t))/Cf + (omega_pll + omega_0)*omega_b*Vq_t               # Non-state variable - valid
    Vq_t_dot    = (omega_b*(Iq_s - Iq_t))/Cf - (omega_pll + omega_0)*omega_b*Vd_t               # Non-state variable - valid
    
    fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
    
    return fX
end
