#Tutorials on PowerSystemsdyanamics.jl
#https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/blob/main/docs/src/tutorials/tutorial_240bus.md

using PowerSimulationsDynamics
using PowerSystemCaseBuilder
using PowerSystems
const PSY = PowerSystems
using Sundials
using Plots
using OrdinaryDiffEq

# %%

# exploring test cases in the library.

using PowerSystemCaseBuilder
show_systems()

# %%

# We remove the checks in this example to avoid large prints
sys = build_system(PSIDSystems, "14 Bus Base Case"; runchecks = false)

# Transform the system's load
for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
end

# %%

using Logging
sim_ida = Simulation(
    ResidualModel,
    sys, #system
    pwd(),
    (0.0, 20.0), #time span
    BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
    console_level = Logging.Info,
)

execute!(sim_ida, IDA(), dtmax = 0.01)

# %%

res_ida = read_results(sim_ida)
v1101_ida = get_voltage_magnitude_series(res_ida, 1101);
plot(v1101_ida)

# %%
sim_rodas = Simulation(
    MassMatrixModel,
    sys, #system
    pwd(),
    (0.0, 20.0), #time span
    BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
    console_level = Logging.Info,
)

res_rodas = read_results(sim_rodas)

# %%


execute!(
    sim_rodas,
    Rodas4(),
    saveat = 0.01,
    atol = 1e-10,
    rtol = 1e-10,
    initializealg = NoInit(),
)

#%%

v1101 = get_voltage_magnitude_series(res_rodas, 1101);
plot(v1101, label = "RODAS4")
plot!(v1101_ida, label = "IDA")
