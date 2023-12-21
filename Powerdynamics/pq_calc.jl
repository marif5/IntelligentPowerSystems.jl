wc = 50;
Vt = 0.9;
theta_t = pi/4;
theta_pll = 10*(pi/180);

Vr = Vt*cos(theta_t);
Vi = Vt*sin(theta_t);

It = (Vt*exp(theta_t) - 1)/(im*0.1);
Ir = abs(It)*cos(angle(It));
Ii = abs(It)*sin(angle(It));

Vdq = [cos(theta_pll) -1*sin(theta_pll); sin(theta_pll) cos(theta_pll)]*[Vr; Vi];
Idq = [cos(theta_pll) -1*sin(theta_pll); sin(theta_pll) cos(theta_pll)]*[Ir; Ii];

p = Vdq[1,:]*Idq[1,:]' + Vdq[2,:]*Idq[2,:]';
q = Vdq[2,:]*Idq[1,:]' - Vdq[1,:]*Idq[2,:]';

p_fil = ()