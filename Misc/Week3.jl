# https://docs.juliahub.com/PowerSimulationsDynamics/T1QyN/0.1.2/Examples/example_data/

#using Pkg
#Pkg.activate(.)
#Pkg.add("PowerSystemDynamics")
#Pkg.add("PowerSystem")

using PowerSimulationsDynamics
using PowerSystems
using Sundials
const PSY = PowerSystems

# %% ==============


omib_sys = System("omib_sys.json")

#Compute the Y_bus after fault
#Collect the branch of the system as:
fault_branch = deepcopy(collect(get_components(Branch, omib_sys))[1])
#Duplicates the impedance of the reactance
fault_branch.x = fault_branch.x * 2
#Obtain the new Ybus of the faulted system
Ybus_fault = Ybus([fault_branch], get_components(Bus, omib_sys))[:, :]

#Construct the perturbation
perturbation_Ybus = ThreePhaseFault(
    1.0, #change will occur at t = 1.0s
    Ybus_fault, #new Ybus
)

#Time span of our simulation
tspan = (0.0, 30.0)

#Define Simulation
sim = Simulation(
    pwd(),
    sys, #system
    tspan, #time span
    perturbation_Ybus, #Type of perturbation
)

#Will print the initial states. It also give the symbols used to describe those states.
print_device_states(sim)
#Will export a dictionary with the initial condition values to explore
x0_init = get_initial_conditions(sim)

#Solve problem
run_simulation!(sim, #simulation structure
                IDA(), #Sundials DAE Solver
                dtmax=0.02); #Arguments: Maximum timestep allowed

                using Plots
angle = get_state_series(sim, ("generator-102-1", :δ))
plot(angle, xlabel="time", ylabel="rotor angle [rad]", label="rotor angle")

volt = get_voltagemag_series(sim, 102)
plot(volt, xlabel="time", ylabel="Voltage [pu]", label="V_2")

d   