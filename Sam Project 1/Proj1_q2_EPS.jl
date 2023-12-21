# ================= Q2 -Newton raphson warmup ==============#
# ==========================================================#
using Plots

# Function definition
function NWR(x1_0::Float64, x2_0::Float64, m::Float64)
        i = 1;
        e1 = 0.01;
        e2 = 0.01;
        er_x1 = zeros(110, 1);
        er_x2 = zeros(110, 1);
        x1 = zeros(110, 1);
        x2 = zeros(110, 1);

        x1[1] = x1_0;
        x2[1] = x2_0;
        er_x1[1] = e1;
        er_x2[1] = e2;

        while abs(er_x1[i]) > 0.005 || abs(er_x2[i]) > 0.005
        
            Fi = [exp(x1[i]*x2[i]) ; cos(x1[i]+x2[i])]
            Ji = [x2[i]*exp(x1[i]*x2[i]) x1[i]*exp(x1[i]*x2[i]) ; -1*sin(x1[i]+x2[i]) -1*sin(x1[i]+x2[i])]

            X = [x1[i] ; x2[i]] - inv(Ji)*Fi
            x1[i+1] = X[1,1]
            x2[i+1] = X[2,1]

           if m == 1.0
                er_x1[i+1] = maximum(abs.(x1[i+1] - x1[i]) ./ x1[i])
                er_x2[i+1] = maximum(abs.(x2[i+1] - x2[i]) ./ x2[i])
           
            elseif m == 2.0
                
                er_x1[i+1] = maximum(Fi[1, 1])
                er_x2[i+1] = maximum(Fi[2, 1])
            end
            
            i = i + 1    
        end

    return x1, x2, er_x1, er_x2, i
end

x1, x2, er_x1, er_x2, i= NWR(1.2,0.5,1.0)
println("Final values:")
println("x1 = ", x1[i])
println("x2 = ", x2[i])

plot([1:i],x1[1:i])
plot!([1:i],x2[1:i])
savefig("2_1.png")

plot([1:i],er_x1[1:i])
plot!([1:i],er_x2[1:i])
savefig("2_2.png")