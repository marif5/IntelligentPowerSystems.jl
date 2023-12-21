using SymPy, JuMP, LinearAlgebra, DataFrames, CSV, Plots, Ipopt

# Including fucntions file
include("UGF_funcs.jl")

# %% MAIN
#====================================================================================#

# Metered values at the terminal/PCC

i = 1
Vabc = zeros(3, 11) 
theta_t = zeros(11) 
for t in 0.0:0.001:0.01 
    Vabc[:, i] = [1.0 * sin(2 * pi * 60 * t);            # Phase A
                  1.0 * sin(2 * pi * 60 * t - 2*pi/3);   # Phase B
                  1.0 * sin(2 * pi * 60 * t + 2*pi/3)]   # Phase C
    theta_t[i] = 2 * pi * 60 * t + 0.25
    local ind = i
    i += 1
end

Vdq0 = parks_transf(Vabc, theta_t)
Vdt = Vdq0[2,1]
Vqt = Vdq0[1,1]


# %% Defining coefficents

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
V_inf = 1
theta_inf = 0
V_t0, theta_t0, Vs0, delta0 = init_inv(p0, q0, V_inf, theta_inf , omega_0, Lf, Cf, line_X)

# %% State vector formation

X = zeros(12);
    X[1] = p0                                   #P_hat      State + (Eq5)
    X[2] = q0                                   #Q_hat      State + (Eq6)
    X[3] = 0                                    #zeta       State + (Eq7)
    X[4] = theta_t0                             #theta_pll          IPOPT
    X[5] = delta0                               #delta              IPOPT
    X[6] = 0                                    #phi_d      State + (Eq14)
    X[7] = 0                                    #gamma_d    State + (Eq16)
    X[8] = Vs0                                  #Vs_abs             IPOPT

        Vd_t = V_t0*cos(theta_t0)                                          
        Vq_t = V_t0*sin(theta_t0)
        Id_t = (Vd_t-V_inf)/line_X                                  # Non-state variable - valid   
        Iq_t = Vq_t/line_X

        omega_pll_0 = Kp_pll*X[3]                                 #  (Eq8) 
        Id_s    = -1*(omega_pll_0 + omega_0)*Vq_t*Cf + Id_t       #  (Eq23)
        Iq_s    =  (omega_pll_0 + omega_0)*Vd_t*Cf + Iq_t         #  (Eq24)

    X[9]  = Id_s                                 #Id_s                    
    X[10] = Iq_s                                 #Iq_s                    
    X[11] = Vd_t                                 #Vd_t       (Eq1) + IPOPT
    X[12] = Vq_t                                 #Vq_t       (Eq1) + IPOPT
    
PLL_coeff = [Kp_pll, Ki_pll, omega_0, omega_DQ, omega_b]
P_droop_Coeff = [omega_c, Kp_i, p0, q0, v0, mp, mq]
IVControl_coeff = [Kp_vc, Ki_vc, Kf_vc, Kp_ic, Ki_ic, Kf_ic, Cf, Lf, Kv_i]


# %% Forward Euler solution
include("UGF_funcs.jl")
times, states, OUTX = inverter_der(X, line_X, V_inf, PLL_coeff, P_droop_Coeff, IVControl_coeff)




#=== Debug space======#




#stating coefficents:
Kp_pll      = PLL_coeff[1]
Ki_pll      = PLL_coeff[2]
omega_0     = PLL_coeff[3]
omega_DQ    = PLL_coeff[4]
omega_b     = PLL_coeff[5]

omega_c     = P_droop_Coeff[1]
Kp_i        = P_droop_Coeff[2]
p0          = P_droop_Coeff[3]
q0          = P_droop_Coeff[4]
v0          = P_droop_Coeff[5]
mp          = P_droop_Coeff[6]
mq          = P_droop_Coeff[7]

Kp_vc       = IVControl_coeff[1]
Ki_vc       = IVControl_coeff[2]
Kf_vc       = IVControl_coeff[3]
Kp_ic       = IVControl_coeff[4]
Ki_ic       = IVControl_coeff[5]
Kf_ic       = IVControl_coeff[6]
Cf          = IVControl_coeff[7]
Lf          = IVControl_coeff[8]
Kv_i        = IVControl_coeff[9]

Xi = copy(X)

# Forward Euler formation    
    h = 0.02    
    times = 0.02:h:0.4
   # for i in 1:3
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
    #Id_t = (Vd_t-V_inf0)/xline                                 # Non-state variable - valid   
    #Iq_t = Vq_t/xline                                          # Non-state variable - valid

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
    Id_s_dot    = (omega_b*(Vd_s - Vd_t))/Lf + (omega_pll + omega_0)*omega_b*Iq_s      # Non-state variable - valid
    Iq_s_dot    = (omega_b*(Vq_s - Vq_t))/Lf - (omega_pll + omega_0)*omega_b*Id_s      # Non-state variable - valid
    Vd_t_dot    = (omega_b*(Id_s - Id_t))/Cf + (omega_pll + omega_0)*omega_b*Vq_t      # Non-state variable - valid
    Vq_t_dot    = (omega_b*(Iq_s - Iq_t))/Cf - (omega_pll + omega_0)*omega_b*Vd_t      # Non-state variable - valid
    
    
    
    fX = [p_hat_dot, q_hat_dot, zeta_dot, theta_pll_dot, delta_dot, phi_d_dot, gamma_d_dot, Vs_abs_dot, Id_s_dot, Iq_s_dot, Vd_t_dot, Vq_t_dot]
    X  = X + h* fX
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
Out_X = DataFrame(P_hat = Xi[1], Q_hat = Xi[2], zeta = Xi[3], theta_t = Xi[4], delta = Xi[5], phi_d = Xi[6], gamma_d = Xi[7], Vs_abs = Xi[8], Ids = Xi[9], Iqs = Xi[10], Vdt = Xi[11], Vqt = Xi[12])
    print(X)
        X  = X + h* fX

    X1 = X2 = X3 = X4 = X5 = X6 = X6 = X7 = X8 = X9 = X10 = X11 = X12 = zeros(length(times))



   
Out_X = DataFrame(P_hat = X1, Q_hat = X2, zeta = X3, theta_t = X4, delta = X5, phi_d = X6, gamma_d = X7, Vs_abs = X8, Ids = X9, Iqs = X10, Vdt = X11, Vqt = X12)

