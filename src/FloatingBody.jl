"""
    FloatingBody

**FloatingBody.jl** is a numerical simulation for floating bodies (6dof). 

See the documentation for usage in your own work.
"""
module FloatingBody

# List of constants 
const γ_GAS = 1.4
const P_ATM = 101325.0
const ρ_AIR = 1.2
const ρ_WATER = 1025.0
const GRAV = 9.80665 # 9.81


# WAMIT SCALING FACTORS
const TRANSLATIONAL = Set([1, 2, 3, 9, 15, 21, 27, 33])
const ROTATIONAL = Set([4, 5, 6])

@inline function get_scaling_factor_Aij(i, j, rho, L)
    k = (i in TRANSLATIONAL && j in TRANSLATIONAL) ? 3 :
        (i in ROTATIONAL && j in ROTATIONAL) ? 5 : 4
    return rho * L^k
end
@inline function get_scaling_factor_Xi(i, rho, g, L)
    m = i in TRANSLATIONAL ? 2 : 3
    return rho * g * (L^m)
end

# import @reexport now to make it available for further imports/exports
using Reexport: @reexport

using MuladdMacro: @muladd
using Random: Xoshiro, MersenneTwister

using StaticArrays: StaticArrays, MVector, MArray, SMatrix,
    @MVector, @SMatrix
@reexport using StaticArrays: SVector, @SVector

# using Rotations: RotXYZ, UnitQuaternion 

# @reexport using Parameters: @consts
using SimpleUnPack: @pack!
@reexport using SimpleUnPack: @unpack

using Static: Static, One, True, False

using HDF5

using LinearAlgebra

using Dates
using Statistics
using Printf

using SciMLBase: SciMLBase
using OrdinaryDiffEqCore: ODEIntegrator


# Fast, no checks version
function step_mooring_solver_fast!(integrator::ODEIntegrator)
    lines_class = integrator.opts.callback.discrete_callbacks[1].affect!.coupling_manager.lines_class
    bc_times = lines_class[1].bc_ext.boundary_condition_itp.times
    end_time_bc = bc_times[end]
    time_remaining = end_time_bc - integrator.t
    SciMLBase.step!(integrator, time_remaining, true)
    return nothing
end
# Validates variables before updating
function step_mooring_solver!(integrator::ODEIntegrator;
    dt::Real=zero(eltype(integrator.sol)))
    lines_class = integrator.opts.callback.discrete_callbacks[1].affect!.coupling_manager.lines_class
    bc_times = lines_class[1].bc_ext.boundary_condition_itp.times
    @inbounds for line in Iterators.drop(lines_class, 1)
        if line.bc_ext.boundary_condition_itp.times != bc_times
            throw(ArgumentError("Not all line BCs are synchronised: line_id = $(line.line_id)"))
        end
    end
    start_time_bc = bc_times[1]
    end_time_bc = bc_times[end]
    time_remaining = end_time_bc - integrator.t
    if iszero(dt)
        dt_final = time_remaining
    else
        dt_final = dt
        if dt_final > time_remaining
            if time_remaining <= 0.0
                @warn "BC window exhausted" current_time=integrator.t bc_window=(start_time_bc,
                    end_time_bc)
                return nothing
            else
                dt_final = time_remaining
                @warn "Timestep adjusted to BC window" requested_dt=dt actual_dt=dt_final
            end
        end
        dt_final = min(dt, time_remaining)
    end
    SciMLBase.step!(integrator, dt_final, true)
    return nothing
end

# Dispatch update of boundary conditions from the ODE integrator.
function set_boundary_conditions!(integrator::ODEIntegrator,
    coupling_id::Int, line_id::Int,
    new_times, new_u)
    external_callback = integrator.opts.callback.discrete_callbacks[1]
    set_boundary_conditions!(external_callback.affect!,
        coupling_id, line_id,
        integrator.t,
        new_times, new_u)
end
function set_boundary_conditions!(integrator::ODEIntegrator, coupling_id::Int, new_times,
    new_u)
    external_callback = integrator.opts.callback.discrete_callbacks[1]
    set_boundary_conditions!(external_callback.affect!,
        coupling_id, 1, # only 1 line assumed
        integrator.t,
        new_times, new_u)
end
@inline function set_boundary_conditions!(external_coupling,
    coupling_id::Int,
    line_id::Int,
    t::Real,
    new_times, # ::AbstractVector,
    new_u) #::AbstractMatrix)
    manager = external_coupling.coupling_manager

    # Update interpolator
    line = manager.lines_class[line_id]
    boundary_condition_itp = line.bc_ext.boundary_condition_itp
    unsafe_set_boundary_conditions!(boundary_condition_itp, t, new_times, new_u)

    return nothing
end
@inline function unsafe_set_boundary_conditions!(interp, t,
    new_times::RealT,
    new_u::AbstractVector{RealT}) where {
    RealT
}
    # Update the time vector (only 2 elements)
    interp.times[1] = interp.times[2]
    interp.times[2] = new_times

    @inbounds for i in 1:size(interp.u, 1)
        interp.u[i, 1] = interp.u[i, 2]
        interp.u[i, 2] = new_u[i]
    end
    return nothing
end
# Get mooring forces
function get_mooring_force(integrator::ODEIntegrator, coupling_id::Int=1)
    external_callback = integrator.opts.callback.discrete_callbacks[1]
    return get_mooring_force(external_callback.affect!, coupling_id)
end
@inline function get_mooring_force(external_coupling,
    coupling_id::Int)
    return external_coupling.coupling_manager.mooring_3dforces
end

# v1, v3, v5, v9, v15, v21, v27, v33 = states[1:8]
# z1, z3, z5, z9, z15, z21, z27, z33 = states[9:16]

# p_rel = states[17:21]
# Omega = states[22:26]
# E_turb = states[27:31]

#=
  Include all top-level source files
=#
include("wave_spectra.jl")
include("utilities.jl")
include("force_radiation.jl")
include("force_excitation.jl")
include("pto.jl")
include("buoy_simul.jl")
include("rhs_funcs.jl")
include("body_rotations.jl")
include("output.jl")

export WaveClimateAtlantic, rhs_f!, rhs_f_static!,
    BiradialTurbine, OWC_piston, FloatingBodySim, SimulationOutputs, WaveSpectra
export nvars, nmodes, nowcs
export compute_output!

end
