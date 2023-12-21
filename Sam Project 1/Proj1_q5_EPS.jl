using LinearAlgebra
using CSV
using Plots
using DataFrames

path_1 = "C:\\Users\\marif\\OneDrive - University of Vermont\\Documents\\PhD UVM\\Research\\VScode\\IntelligentPowerSystems.jl\\IntelligentPowerSystems.jl\\BranchData.csv"
path_2 = "C:\\Users\\marif\\OneDrive - University of Vermont\\Documents\\PhD UVM\\Research\\VScode\\IntelligentPowerSystems.jl\\IntelligentPowerSystems.jl\\BusData.csv"
linedata = CSV.File(path_1) |> DataFrame
BusData = CSV.File(path_2) |> DataFrame

linedata.Z = linedata.R + im * linedata.X
linedata.Y = 1.0 ./ linedata.Z

bus = size(BusData, 1)
line = size(linedata, 1)
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
Y
spy(abs.(Y))


# Define initial values for voltage magnitudes and angles and eliminating slack bus parameters
#V1, ang1, V2 = 1.0, 0.0, 1.0

X = hcat(BusData.Bus_delta, BusData.Bus_mag) 
X = X[setdiff(1:size(X, 1), 21), :] 
X = DataFrame(X, [:V_delta, :V_mag])

P = zeros(Complex, bus+1, 1)

# Power equations
 
for k in 1:bus, 
    for n in 1:bus
    P[k+1,1] =  P[k,1] + abs(Y[k, n]) *  X.V_mag[n]  * cos(X.V_delta[k] - X.V_delta[n] - angle(Y[k, n]))
    end
    P[k+1,1] = X.V_mag[k] * P[k+1,1]
end

P1 = V1 * V1 * abs(Y[1, 1]) * cos(ang1 - ang1 - angle(Y[1, 1])) + V1 * V2 * abs(Y[1, 2]) * cos(ang1 - ang2 - angle(Y[1, 2])) + V1 * V3 * abs(Y[1, 3]) * cos(ang1 - ang3 - angle(Y[1, 3]))
P2 = V2 * V2 * abs(Y[2, 2]) * cos(ang2 - ang2 - angle(Y[2, 2])) + V1 * V2 * abs(Y[2, 1]) * cos(ang2 - ang1 - angle(Y[1, 2])) + V2 * V3 * abs(Y[2, 3]) * cos(ang2 - ang3 - angle(Y[2, 3]))
P3 = V3 * V3 * abs(Y[3, 3]) * cos(ang3 - ang3 - angle(Y[3, 3])) + V3 * V2 * abs(Y[3, 2]) * cos(ang3 - ang2 - angle(Y[3, 2])) + V1 * V3 * abs(Y[3, 1]) * cos(ang3 - ang1 - angle(Y[3, 1]))

Q1 = V1 * V1 * abs(Y[1, 1]) * sin(ang1 - ang1 - angle(Y[1, 1])) + V1 * V2 * abs(Y[1, 2]) * sin(ang1 - ang2 - angle(Y[1, 2])) + V1 * V3 * abs(Y[1, 3]) * sin(ang1 - ang3 - angle(Y[1, 3]))
Q2 = V2 * V2 * abs(Y[2, 2]) * sin(ang2 - ang2 - angle(Y[2, 2])) + V2 * V1 * abs(Y[2, 1]) * sin(ang2 - ang1 - angle(Y[2, 1])) + V2 * V3 * abs(Y[2, 3]) * sin(ang2 - ang3 - angle(Y[2, 3]))
Q3 = V3 * V3 * abs(Y[3, 3]) * sin(ang3 - ang3 - angle(Y[3, 3])) + V3 * V2 * abs(Y[3, 2]) * sin(ang3 - ang2 - angle(Y[3, 2])) + V3 * V1 * abs(Y[3, 1]) * sin(ang3 - ang1 - angle(Y[3, 1]))

