# Import necessary libraries
using LinearAlgebra
using DataFrames
using SymPy

# Define line data and bus data using DataFrames
linedata = DataFrame(fbus=[1, 1, 2], tbus=[2, 3, 3], R=[0.05, 0.05, 0.05], X=[0.1, 0.1, 0.1])
line = size(linedata, 1)

busdata = DataFrame(bus=[1, 2, 3], Pd=[1.9, 0.8, -1.8], Qd=[0, 0, 0])

# Calculate Z and Y parameters from R and X
linedata.Z = linedata.R + im * linedata.X
linedata.Y = 1.0 ./ linedata.Z

# Get the number of buses
bus = size(busdata, 1)

# Create a matrix Y to represent the admittance matrix
Y = zeros(Complex, bus, bus)

# Fill in the non-diagonal elements
for l in 1:line
    Y[linedata[l, :fbus], linedata[l, :tbus]] = -linedata[l, :Y]
    Y[linedata[l, :tbus], linedata[l, :fbus]] = -linedata[l, :Y]
end

# Fill in the diagonal elements
for b in 1:bus, l in 1:line
    if (b == linedata[l, :fbus]) || (b == linedata[l, :tbus])
        Y[b, b] += linedata[l, :Y]
    end
end

# Define initial values for voltage magnitudes and angles
V1, ang1, V2 = 1.0, 0.0, 1.0

# Symbolic variables for voltage magnitudes and angles of buses 2 and 3
V3, ang2, ang3 = symbols("V3, ang2, ang3", real=true)

# Power equations
P1 = V1 * V1 * abs(Y[1, 1]) * cos(ang1 - ang1 - angle(Y[1, 1])) + V1 * V2 * abs(Y[1, 2]) * cos(ang1 - ang2 - angle(Y[1, 2])) + V1 * V3 * abs(Y[1, 3]) * cos(ang1 - ang3 - angle(Y[1, 3]))
P2 = V2 * V2 * abs(Y[2, 2]) * cos(ang2 - ang2 - angle(Y[2, 2])) + V1 * V2 * abs(Y[2, 1]) * cos(ang2 - ang1 - angle(Y[1, 2])) + V2 * V3 * abs(Y[2, 3]) * cos(ang2 - ang3 - angle(Y[2, 3]))
P3 = V3 * V3 * abs(Y[3, 3]) * cos(ang3 - ang3 - angle(Y[3, 3])) + V3 * V2 * abs(Y[3, 2]) * cos(ang3 - ang2 - angle(Y[3, 2])) + V1 * V3 * abs(Y[3, 1]) * cos(ang3 - ang1 - angle(Y[3, 1]))

Q1 = V1 * V1 * abs(Y[1, 1]) * sin(ang1 - ang1 - angle(Y[1, 1])) + V1 * V2 * abs(Y[1, 2]) * sin(ang1 - ang2 - angle(Y[1, 2])) + V1 * V3 * abs(Y[1, 3]) * sin(ang1 - ang3 - angle(Y[1, 3]))
Q2 = V2 * V2 * abs(Y[2, 2]) * sin(ang2 - ang2 - angle(Y[2, 2])) + V2 * V1 * abs(Y[2, 1]) * sin(ang2 - ang1 - angle(Y[2, 1])) + V2 * V3 * abs(Y[2, 3]) * sin(ang2 - ang3 - angle(Y[2, 3]))
Q3 = V3 * V3 * abs(Y[3, 3]) * sin(ang3 - ang3 - angle(Y[3, 3])) + V3 * V2 * abs(Y[3, 2]) * sin(ang3 - ang2 - angle(Y[3, 2])) + V3 * V1 * abs(Y[3, 1]) * sin(ang3 - ang1 - angle(Y[3, 1]))

# Jacobian matrix
J11 = [diff(P2, ang2) diff(P2, ang3); diff(P3, ang2) diff(P3, ang3)]
J12 = [diff(P2, V3); diff(P3, V3)]
J21 = [diff(Q3, ang2) diff(Q3, ang3)]
J22 = [diff(Q3, V3)]
J = [J11 J12; J21 J22]

# Newton-Raphson iteration loop
function NR_loop()
    tovars = [ang2, ang3, V3]
    toval = [0.0, 0.0, 1.0]
    iter = 0
    gan = 10
    while gan > 1e-3 && iter < 100
        Pcalc = [P2.subs(zip(tovars, toval)) P3.subs(zip(tovars, toval))]
        Qcalc = Q3.subs(zip(tovars, toval))
        Pd = busdata.Pd[2:end]
        y_delta = [Pd[1] Pd[2] 0] - [Pcalc Qcalc]
        J_calc = J.subs(zip(tovars, toval))
        x_delta = inv(J_calc) * (y_delta)'
        toval = toval + x_delta
        iter = iter + 1
        gan = norm(y_delta)
    end
    return toval, gan, iter
end

# Execute the Newton-Raphson loop
results = NR_loop()
