using FloatingBody, OrdinaryDiffEq, OrdinaryDiffEqSSPRK

dtype = Float64
modes = [1, 3, 5, 9, 15, 21, 27, 33]
wave_spec = WaveSpectra{dtype}(WaveClimateAtlantic.HS_TABLE[1], WaveClimateAtlantic.TP_TABLE[1], n_waves=400, ramp=20.0)

# Instatiate turbine
turb1 = BiradialTurbine{dtype}(1.2; pfb=3.0, Psi_max=10.0, Psi_bep=0.3570984495110234, I_turb_ref=5.01, P_rated=20000.0)

# Instatiate PTO
owc9 = OWC_piston{dtype}(9, chamber_height=15.0, chamber_area=304.34, n_turbines=5, turbine=turb1)
owc15 = OWC_piston{dtype}(15, chamber_height=15.0, chamber_area=481.25, n_turbines=5, turbine=turb1)
owc21 = OWC_piston{dtype}(21, chamber_height=15.0, chamber_area=481.25, n_turbines=5, turbine=turb1)
owc27 = OWC_piston{dtype}(27, chamber_height=15.0, chamber_area=481.25, n_turbines=5, turbine=turb1)
owc33 = OWC_piston{dtype}(33, chamber_height=15.0, chamber_area=481.25, n_turbines=5, turbine=turb1)
vec_PTOs = [owc9, owc15, owc21, owc27, owc33]

# Instantiate OCTAPLAT
fbody_siml = FloatingBodySim{dtype}("src/OCTAPLAT.h5", "src/radiation_models.h5", modes, vec_PTOs, wave_spec)



# Initial conditions
Ω_0 = 70.0
state_var_n = 3 + 5
# Turbines start and end at:
i_turb = state_var_n * 2 + 1 # 17

#!!
u0 = zeros(dtype, nvars(fbody_siml))
u0[i_turb+1] = Ω_0
u0[i_turb+1+3] = Ω_0
u0[i_turb+1+3+3] = Ω_0
u0[i_turb+1+3+3+3] = Ω_0
u0[i_turb+1+3+3+3+3] = Ω_0

# du0 = zeros(6)
dt = 0.005
tspan = (0.0, 1000)
prob = ODEProblem(rhs_f!, u0, tspan, fbody_siml)
sol = solve(prob, SSPRK43(), dt=dt)

using Plots
u_sol = sol.u
t_sol = sol.t
n_timesteps = length(t_sol)
output = SimulationOutputs{dtype}(n_timesteps, nmodes(fbody_siml), nowcs(fbody_siml))
compute_output!(output, sol, fbody_siml)
plot(output.time, output.velocity[:,2], title = "velocity heave")

# v3 = map(x -> x[2], sol.u)
# plot(v3)

# integ = init(prob, SSPRK33(), dt=dt)
# step!(integ)