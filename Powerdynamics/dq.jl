using Plots

w = 2 * pi * 60
Vrms = (1/sqrt(2)) * (12.47 * 1000)
theta_t = 10*pi/180

t_values = 0:0.0001:(1/60)

Va = Complex{Float64}[]
Vb = Complex{Float64}[]
Vc = Complex{Float64}[]

    gif_Vabc = @animate for t in t_values
        Va_theta = Vrms * exp(im * (w * t + theta_t))
        Vb_theta = Vrms * exp(im * (w * t + theta_t - 2*pi/3))
        Vc_theta = Vrms * exp(im * (w * t + theta_t + 2*pi/3))
        
        push!(Va, Va_theta)
        push!(Vb, Vb_theta)
        push!(Vc, Vc_theta)
        
        plot([0, real(Va_theta)], [0, imag(Va_theta)], label="Va", linecolor=:blue, legend=true, xlims=(-15000, 15000), ylims=(-10000, 10000))
        plot!([0, real(Vb_theta)], [0, imag(Vb_theta)], label="Vb", linecolor=:blue)
        plot!([0, real(Vc_theta)], [0, imag(Vc_theta)], label="Vc", linecolor=:blue)
        xlabel!("Real Part")
        ylabel!("Imaginary Part")
        title!("ABC Complex Voltage")
    end

gif(gif_Vabc, "Vabc.gif", fps = 60)
V_abc = [ Va';  Vb';  Vc'];

A_dq0 = [cos(theta_t) sin(theta_t) 0; -sin(theta_t) cos(theta_t) 0; 0 0 1];
A_c = (sqrt(2) / 3)*[1.0 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2; 1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
A_p = A_dq0 * A_c;

V_dq0 = A_p * V_abc


    gif_Vdq0 = @animate for i in 1:length(t_values)
        Vd=V_dq0[1,i]
        Vq=V_dq0[2,i]
        V0=V_dq0[3,i]
            plot([0, real(Vd)], [0, imag(Vd)], label="Vd", linecolor=:red, legend=true, xlims=(-15000, 15000), ylims=(-10000, 10000))
            plot!([0, real(Vq)], [0, imag(Vq)], label="Vq", linecolor=:red)
            plot!([0, real(V0)], [0, imag(V0)], label="V0", linecolor=:red)
            xlabel!("Real Part")
            ylabel!("Imaginary Part")
           title!("DQ0 Complex Voltage")
    end
    
gif(gif_Vdq0, "Vdq0.gif", fps = 60)
plot([real(V_dq0[1,:], imag(V_dq0[1,:])

#Playing with A_p (Some other time)

#Power Calculations and filtering
#Assumiptions: Single bus with slack case

Vd = V_dq0[1,:]
Vq = V_dq0[2,:]
