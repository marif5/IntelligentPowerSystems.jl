#Assuming 12.47kV as the predominant distribution voltage in VT

vt_a = ((12.47*1000)*cos(5*pi/180)+((12.47*1000)*sin(5*pi/180))*im ) / (1.732);
vt_b = vt_a*cis(4π/3);
vt_c = vt_a*cis(2π/3); 
vt_base = 12.47*1000/1.732;

theta_pll = 10*pi/180;          #assumption ??

#Parkstransformation
VABC = [vt_a vt_b vt_c]'
P = 2/3*[sin(theta_pll) sin]
VDQ0 = VABC;


wb = 2*pi*60;
#global (DQ) reference parameters

vt_D = 1;
vt_Q = 1;
it_D = 1;
it_Q = 1;

vt_DQ = vt_D + vt_Q*im;
it_DQ = it_D + it_Q*im;
theta_t = angle(vt_DQ)     #theta_t is the phase of terminal voltage in DQ refernce

#local (dq) inverter reference parameters

vt_dq = 1 + 1*im;
it_dq = 0.1 + 0.1*im;
wc = 1;                 #wc is the lowspass filter cut-off freq in per unit
F = [0.9; 0.09];        #F is the transformation matrix for low pass filer

vt_d = real(it_dq);
vt_q = imag(it_dq);
it_d = real(it_dq);
it_q = imag(it_dq);

# S is the dq-framed computed (instantaneous) power for the inverter 'i'
# S_hat is the lowpass filtered version of p and q

#S = [vt_d vt_q; vt_q -vt_d]*([it_d it_q]');
#S_hat = F.*S;
#p_hat = wc*(S[1] - S_hat[1]);
#q_hat = wc*(S[2] - S_hat[2]);


