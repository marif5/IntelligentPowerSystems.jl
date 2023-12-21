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
    Y[linedata[l, :fbus], linedata[l, :tbus]] = (-BranchData[l,:Y])./(conj(BranchData[l,:ont]))
    Y[linedata[l, :tbus], linedata[l, :fbus]] = -BranchData[l,:Y]./(BranchData[l,:ont])

end

# Diagonal elements
for b in 1:Buses, l in 1:Lines
    if (b==BranchData[l,:fbus])
        Y[b,b] += (BranchData[l,:Y] + 0.5*im*BranchData[l,:b] + BusData[b,:Bsh])./(abs.(BranchData[l,:ont]))^2
        (b==BranchData[l,:tbus])
    end
    if (b==BranchData[l,:tbus])
        Y[b,b] += (BranchData[l,:Y] + 0.5*im*BranchData[l,:b] + BusData[b,:Bsh])
    end
end

#================= Q3-b: Phase shifting transformers =======================#

ph_shft_b71_b73 = angle(0.9962+0.0872*im)

#================= Q3-c: Spy command (Sparsity and Symmetry)================#

spy(abs.(Y))
savefig("3_a1.png")
spy(abs.(Y'))
savefig("3_a2.png")

#=================Q3-d Creating sub-matrix out of Nth order ===============#

N = 5
last_rows = size(Y, 1) - (N-1):size(Y, 1)
last_columns = size(Y, 2) - (N-1):size(Y, 2)
sub_Y = Y[last_rows, last_columns]

spy(abs.(sub_Y))
savefig("3_d1.png")

#=================Q3-e Weighted Admittance ================#

from = linedata.fbus
to = linedata.tbus
E = zeros(bus,line)

for f in 1:bus,t in 1:line
    if (from[f] == f) && (from[t] == t)
        E[f,t] += +1
        E[t,f] += -1
    end
end

#Incpmplete - Faulty results

#===========================================================#
#=================Q7- DC Power Flow=========================#

Pdiff_new = zeros(Float64, bus)

Pdiff = (BusData.Bus_Pg)./100 - (BusData.Bus_Pd)./100
Pdiff_new = Pdiff[1:end .∉ [21]]
B = imag.(Y)
B = B[1:end .∉ [21], 1:end .∉ [21]]
ang = inv(B)*Pdiff_new
plot(ang)
savefig("7_a.png")