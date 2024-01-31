# Including fucntions file
using JuMP, Ipopt, DifferentialEquations, LinearAlgebra, Plots

include("../Powerdynamics/UGF_funcs.jl")

# %%  === Initialization

line_X      = 0.2;
omega_c     = 50;
Kp_pll      = 0.2;
Ki_pll      = 5.0;
Kp_i        = 0.3;
omega_0     = 1.0;
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
v_t0        = 1
theta_inf   = 0
V_inf, theta_t0, Vs0, delta0 = init_inv(p0, q0, v_t0, theta_inf , omega_0, Lf, Cf, line_X)

# State vector formation
X    = zeros(12);
X[1] = p0                                   #P_hat      State + (Eq5)
X[2] = q0                                   #Q_hat      State + (Eq6)
X[3] = 0                                    #zeta       State + (Eq7)
X[4] = theta_t0                             #theta_pll          IPOPT
X[5] = delta0                               #delta              IPOPT
X[6] = 0                                    #phi_d      State + (Eq14)
X[7] = 0                                    #gamma_d    State + (Eq16)
X[8] = Vs0                                  #Vs_abs             IPOPT

VD_t = v_t0*cos(theta_t0)                                          
VQ_t = v_t0*sin(theta_t0)

Vd_t = v_t0                                        
Vq_t = 0

IQ_t = (V_inf - VD_t)/line_X                              # Non-state variable - valid   
ID_t = VQ_t/line_X                                        # changes

R    = [cos(theta_t0) -sin(theta_t0); sin(theta_t0) cos(theta_t0)]
Idq  = R\[ID_t; IQ_t]
Id_t = Idq[1]
Iq_t = Idq[2]

Vdq  = R\[VD_t; VQ_t]
omega_pll_0 = Kp_pll*X[3]                                 #  (Eq8) 

Id_s    =  Id_t - (omega_pll_0 + omega_0)*Vq_t*Cf         #  (Eq23)
Iq_s    =  Iq_t + (omega_pll_0 + omega_0)*Vd_t*Cf         #  (Eq24)

X[9]  = Id_s                                 #Id_s                    
X[10] = Iq_s                                 #Iq_s                    
X[11] = Vd_t                                 #Vd_t       (Eq1) + IPOPT
X[12] = Vq_t                                 #Vq_t       (Eq1) + IPOPT

PLL_coeff       = [Kp_pll, Ki_pll, omega_0, omega_DQ, omega_b]
P_droop_Coeff   = [omega_c, Kp_i, p0, q0, v0, mp, mq]
IVControl_coeff = [Kp_vc, Ki_vc, Kf_vc, Kp_ic, Ki_ic, Kf_ic, Cf, Lf, Kv_i]

tspan = (0.0, 1.0)

# %% ==

function ODE_inverter(dX, X, p, tspan)
    PLL_coeff       = p[1:5]
    P_droop_Coeff   = p[6:12]
    IVControl_coeff = p[13:21]
   
    Kp_pll      = copy(PLL_coeff[1])
    Ki_pll      = copy(PLL_coeff[2])
    omega_0     = copy(PLL_coeff[3])
    omega_DQ    = copy(PLL_coeff[4])
    omega_b     = copy(PLL_coeff[5])

    omega_c     = copy(P_droop_Coeff[1])
    Kp_i        = copy(P_droop_Coeff[2])
    p0          = copy(P_droop_Coeff[3])
    q0          = copy(P_droop_Coeff[4])
    v0          = copy(P_droop_Coeff[5])
    mp          = copy(P_droop_Coeff[6])
    mq          = copy(P_droop_Coeff[7])

    Kp_vc       = copy(IVControl_coeff[1])
    Ki_vc       = copy(IVControl_coeff[2])
    Kf_vc       = copy(IVControl_coeff[3])
    Kp_ic       = copy(IVControl_coeff[4])
    Ki_ic       = copy(IVControl_coeff[5])
    Kf_ic       = copy(IVControl_coeff[6])
    Cf          = copy(IVControl_coeff[7])
    Lf          = copy(IVControl_coeff[8])
    Kv_i        = copy(IVControl_coeff[9])

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

    R           = [cos(theta_pll) -sin(theta_pll); sin(theta_pll) cos(theta_pll)]
    vdq         = [Vd_t; Vq_t]
    VDQ         = R*vdq                                # DQ to dq (Global to Local)
    VD_t        = VDQ[1]
    VQ_t        = VDQ[2]

    IQ_t        = (V_inf - VD_t)/line_X                # Power balance equation   
    ID_t        = VQ_t/line_X                          # Power balance equation

    Idq         = R\[ID_t; IQ_t]                       # dq to DQ (Local to Global)
    Id_t        = Idq[1]
    Iq_t        = Idq[2]

    p           = Vd_t*Id_t + Vq_t*Iq_t                         
    q           = Vq_t*Id_t - Vd_t*Iq_t                         
    dX[1]       = omega_c*(p-p_hat)                                 # p_hat_dot
    dX[2]       = omega_c*(q-q_hat)                                 # q_hat_dot                        
    theta_t     = atan(VQ_t/VD_t)                               
    
    dX[3]       = theta_t - theta_pll                               # zeta_dot
    omega_pll   = Kp_pll*(theta_t-theta_pll) + Ki_pll*zeta   
    dX[4]       = (omega_pll + omega_0 - omega_DQ)*omega_b           # theta_pll_do

    p_star      = p0 - mp*omega_pll                             
    dX[5]       = Kp_i*(p_star - p_hat)                             # delta_dot
    v_star      = v0 - mq*(q_hat - q0)                          

    dX[6]       = v_star - Vd_t                                     # phi_d_dot
    Id_s_hat    = Kp_vc*(v_star - Vd_t) + Ki_vc*phi_d  + Kf_vc*Id_t - (omega_pll + omega_0)*Cf*Vq_t           # Non-state variable - valid

    dX[7]       = Id_s_hat - Id_s                                     # gamma_d_dot
    Vd_s_hat    = Kp_ic*(Id_s_hat - Id_s  ) + Ki_ic*gamma_d  - Kf_ic*Vd_t - (omega_pll + omega_0)*Lf*Iq_s       # Non-state variable - valid

    Vt_mag      = sqrt((Vd_t)^2 + (Vq_t)^2)
    dX[8]       = Kv_i*(v_star - Vt_mag)                             # Vs_abs_dot
    
    Vd_s        =   Vs_abs*cos(delta - theta_pll)
    Vq_s        =   Vs_abs*sin(delta - theta_pll)

    dX[9]       = (omega_b*(Vd_s - Vd_t))/Lf + (omega_pll + omega_0)*omega_b*Iq_s             # Id_s_dot
    dX[10]      = (omega_b*(Vq_s - Vq_t))/Lf - (omega_pll + omega_0)*omega_b*Id_s             # Iq_s_dot    
    dX[11]      = (omega_b*(Id_s - Id_t))/Cf + (omega_pll + omega_0)*omega_b*Vq_t             # Vd_t_dot         
    dX[12]      = (omega_b*(Iq_s - Iq_t))/Cf - (omega_pll + omega_0)*omega_b*Vd_t             # Vq_t_dot    
    
 end
 
# DX = ODE_inverter(dX, X, tspan, p)
 
 # %% Running ODE Solver

 # ODE_inverter(dX, X, tspan, p)
p = [PLL_coeff; P_droop_Coeff; IVControl_coeff]
# p = Dict("PLL_coeff"         => PLL_coeff,
#         "P_droop_Coeff"      => P_droop_Coeff,    BADDDDDD!! 
#         "IVControl_coeff"    => IVControl_coeff)
# 
# PLL_coeff       = p["PLL_coeff"]
# P_droop_Coeff   = p["P_droop_Coeff"]
# IVControl_coeff = p["IVControl_coeff"]
X0 = copy(X)

prob_ = ODEProblem(ODE_inverter, X0, tspan, p)
sol = solve(prob_, Tsit5())

plot(sol, title = "ODE_inverter", xlabel = "Time", ylabel = "X")

# %%% ==

# %% 
#== Simple Pendulum Problem == #

using DifferentialEquations, Plots

#Constants
const g = 9.81
L = 1.0

#Initial Conditions
u₀ = [0, π / 2]
tspan = (0.0, 6.3)

#Define the problem
function simplependulum(du, u, p, t)
    b = p[1:2]
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)*b[1]
end

#Pass to solvers
p    = [3;4;5]
prob = ODEProblem(simplependulum, u₀, tspan, p)
sol  = solve(prob, Tsit5())

#Plot
plot(sol, linewidth = 2, title = "Simple Pendulum Problem", xaxis = "Time",
    yaxis = "Height", label = ["\\theta" "d\\theta"])

    ==#