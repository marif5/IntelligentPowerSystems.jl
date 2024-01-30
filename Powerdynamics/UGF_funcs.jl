# Fucntion definitions
#====================================================================================#
# DQ Transformation function

# not in use. To be updated and encorporated in rotation (line 98-110)

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

function inverter_dynamics(X, xline, V_inf0, PLL_Coefficient, P_filter_droop_Coefficient, IV_Controller_Coefficient)
    #stating coefficents:
    Kp_pll      = copy(PLL_Coefficient[1])
    Ki_pll      = copy(PLL_Coefficient[2])
    omega_0     = copy(PLL_Coefficient[3])
    omega_DQ    = copy(PLL_Coefficient[4])
    omega_b     = copy(PLL_Coefficient[5])

    omega_c     = copy(P_filter_droop_Coefficient[1])
    Kp_i        = copy(P_filter_droop_Coefficient[2])
    p0          = copy(P_filter_droop_Coefficient[3])
    q0          = copy(P_filter_droop_Coefficient[4])
    v0          = copy(P_filter_droop_Coefficient[5])
    mp          = copy(P_filter_droop_Coefficient[6])
    mq          = copy(P_filter_droop_Coefficient[7])

    Kp_vc       = copy(IV_Controller_Coefficient[1])
    Ki_vc       = copy(IV_Controller_Coefficient[2])
    Kf_vc       = copy(IV_Controller_Coefficient[3])
    Kp_ic       = copy(IV_Controller_Coefficient[4])
    Ki_ic       = copy(IV_Controller_Coefficient[5])
    Kf_ic       = copy(IV_Controller_Coefficient[6])
    Cf          = copy(IV_Controller_Coefficient[7])
    Lf          = copy(IV_Controller_Coefficient[8])
    Kv_i        = copy(IV_Controller_Coefficient[9])

    #stating states:
    p_hat       = copy(X[1])      #powerfilter
    q_hat       = copy(X[2])      #powerfilter
    zeta        = copy(X[3])      #PLL
    theta_pll   = copy(X[4])      #PLL
    delta       = copy(X[5])      #Active Power Control
    phi_d       = copy(X[6])      #IV Contollers
    gamma_d     = copy(X[7])      #IV Contollers
    Vs_abs      = copy(X[8])      #IV Contollers
    Id_s        = copy(X[9])      #IV Contollers - ??
    Iq_s        = copy(X[10])     #IV Contollers - ??
    Vd_t        = copy(X[11])     #IV Contollers - Non Zero
    Vq_t        = copy(X[12])     #IV Contollers - Non Zero

# Explanation for local dq currents:

# Using the local dq axis terminal voltage of the inverter, first the global
# DQ axis terminal voltage are determined by through the rotation matrix R.
# The angle used in 'R' is theta_PLL because its the difference between dq and DQ.
# D and Q axis currents are determined through power equations which are converted
# back to the local dq currents for inverter state equation inputs.

    R   = [cos(theta_pll) -sin(theta_pll); sin(theta_pll) cos(theta_pll)]
    vdq = [Vd_t; Vq_t]
    VDQ = R*vdq                                  # DQ to dq (Global to Local)
    VD_t = VDQ[1]
    VQ_t = VDQ[2]

    IQ_t = (V_inf0 - VD_t)/line_X                # Power balance equation   
    ID_t = VQ_t/line_X                           # Power balance equation

    Idq  = R\[ID_t; IQ_t]                        # dq to DQ (Local to Global)
    Id_t = Idq[1]
    Iq_t = Idq[2]

                            #=
                            theta_t - theta_pll = atan(Vq_t/Vd_t) == theta_t

                            v_t0 = Vd_t
                            VD_t = v_t0*cos(theta_t)                                          
                            VQ_t = v_t0*sin(theta_t)

                            IQ_t = (V_inf0 - VD_t)/line_X                   
                            ID_t = VQ_t/line_X                           

                            R    = [cos(theta_t) -sin(theta_t); sin(theta_t) cos(theta_t)]
                            Idq  = R\[ID_t; IQ_t]
                            Id_t = Idq[1]
                            Iq_t = Idq[2]

                            Vdq  = R\[VD_t; VQ_t]
                            omega_pll_0 = Kp_pll*X[3]                                  

                            Id_s    =  Id_t - (omega_pll_0 + omega_0)*Vq_t*Cf         
                            Iq_s    =  Iq_t + (omega_pll_0 + omega_0)*Vd_t*Cf        

                            =#

    p           = Vd_t*Id_t + Vq_t*Iq_t                         # Non-state variable - valid
    q           = Vq_t*Id_t - Vd_t*Iq_t                         # Non-state variable - valid
    p_hat_dot   = omega_c*(p-p_hat)                             # substate variable - valid
    q_hat_dot   = omega_c*(q-q_hat)                             # substate variable - valid
    theta_t     = atan(VQ_t/VD_t)                               # Non-state variable - valid
    
    zeta_dot       = theta_t - theta_pll                        # Non-state variable - valid
    omega_pll      = Kp_pll*(theta_t-theta_pll) + Ki_pll*zeta   # Non-state variable - valid
    theta_pll_dot  = (omega_pll + omega_0 - omega_DQ)*omega_b   # substate variable - valid

    p_star      = p0 - mp*omega_pll                             # Non-state variable - valid
    delta_dot   = Kp_i*(p_star - p_hat)                         # Non-state variable - valid
    v_star      = v0 - mq*(q_hat - q0)                          # Non-state variable - valid

    phi_d_dot   = v_star - Vd_t                                 # substate variable - valid
    Id_s_hat    = Kp_vc*(v_star - Vd_t) + Ki_vc*phi_d  + Kf_vc*Id_t - (omega_pll + omega_0)*Cf*Vq_t           # Non-state variable - valid

    gamma_d_dot = Id_s_hat - Id_s                               # substate variable - valid
    Vd_s_hat    = Kp_ic*(Id_s_hat - Id_s  ) + Ki_ic*gamma_d  - Kf_ic*Vd_t - (omega_pll + omega_0)*Lf*Iq_s       # Non-state variable - valid

    Vt_mag      = sqrt((Vd_t)^2 + (Vq_t)^2)
    Vs_abs_dot  = Kv_i*(v_star - Vt_mag)                        # Non-state variable - valid
    
    Vd_s    =   Vs_abs*cos(delta - theta_pll)
    Vq_s    =   Vs_abs*sin(delta - theta_pll)

    Id_s_dot    = (omega_b*(Vd_s - Vd_t))/Lf + (omega_pll + omega_0)*omega_b*Iq_s               # Non-state variable - valid
    Iq_s_dot    = (omega_b*(Vq_s - Vq_t))/Lf - (omega_pll + omega_0)*omega_b*Id_s               # Non-state variable - valid
    Vd_t_dot    = (omega_b*(Id_s - Id_t))/Cf + (omega_pll + omega_0)*omega_b*Vq_t               # Non-state variable - valid
    Vq_t_dot    = (omega_b*(Iq_s - Iq_t))/Cf - (omega_pll + omega_0)*omega_b*Vd_t               # Non-state variable - valid
    
    fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
    
    return fX
end





#====================== Archived Fucntions ===================================#

function parkst_abc2dq0(Vabc, theta_t)
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

