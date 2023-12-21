using SymPy, JuMP, LinearAlgebra, DataFrames, CSV, Plots, Ipopt

# Including fucntions file
include("UGF_funcs.jl")

# %% Defining coefficents
include("UGF_funcs.jl")

line_X      = 0.2;
omega_c     = 50;
Kp_pll      = 0.2;
Ki_pll      = 5.0;
Kp_i        = 0.3;
omega_0     = 1;
omega_DQ    = 1.0;
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
V_inf       = 1
theta_inf   = 0
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
Iq_t = (V_inf - Vd_t)/line_X                                  # Non-state variable - valid   
Id_t = Vq_t/line_X

Vd_t*Id_t + Iq_t*Vq_t
Vq_t*Id_t - Vd_t*Iq_t



# %% =======================

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

#  test initialization
include("UGF_funcs.jl")

x_dot = inverter_dynamics(X, line_X, V_inf, PLL_coeff, P_droop_Coeff, IVControl_coeff)
