@muladd begin # Enable FMA
#! format: noindent

function rhs_f_static!(du, u, buoy_sim::FloatingBodySim{NMODES}, t) where NMODES
    @unpack modes = buoy_sim
    @unpack Mg_Inv, C, n_lines = buoy_sim

    # Reset du
    fill!(du, zero(eltype(du)))

    v_vector = get_buoy_vel(u, buoy_sim) # SVector{length(modes)}(@views u[1:length(modes)])
    z_vector = get_buoy_pos(u, buoy_sim)

    f_restore = -(C * z_vector)
    f_pto = calc_pto_force!(du, u, t, buoy_sim)
    f_moor = calc_mooring_force!(du, u, t, buoy_sim)

    du1a8 = Mg_Inv * (f_restore + f_pto + (f_moor .* n_lines))

    @inbounds for (i, j) in enumerate(velrange(buoy_sim))
        du[j] = du1a8[i]
    end

    @inbounds for (i, j) in enumerate(posrange(buoy_sim))
        du[j] = v_vector[i]
    end
    return nothing
end
function rhs_f!(du, u, buoy_sim::FloatingBodySim{NMODES}, t) where NMODES
    @unpack modes = buoy_sim
    @unpack Mg, Mg_Inv, C, n_lines, moor_ramp = buoy_sim

    # Reset du
    fill!(du, zero(eltype(du)))

    v_vector = get_buoy_vel(u, buoy_sim) # SVector{length(modes)}(@views u[1:length(modes)])
    z_vector = get_buoy_pos(u, buoy_sim)

    f_restore = -(C * z_vector)
    # f_exc = calc_excitation_force!(du, u, t, buoy_sim)
    # f_rad = calc_radiation_force!(du, u, t, buoy_sim)
    f_pto = calc_pto_force!(du, u, t, buoy_sim)
    f_moor = calc_mooring_force!(du, u, t, buoy_sim)

    t_ramp = moor_ramp
    if t > t_ramp # ramp dynamic simul
        f_rad = calc_radiation_force!(du, u, t-t_ramp, buoy_sim)
        f_exc = calc_excitation_force!(du, u, t-t_ramp, buoy_sim)

        f_rad_ = f_rad .* SVector(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

        du1a8 = Mg_Inv * (f_restore + f_exc + f_rad_ + f_pto + (f_moor .* n_lines))
    else
        du1a8 = Mg_Inv * (f_restore + f_pto + (f_moor .* n_lines))
    end
    # du1a8 = Mg_Inv * (f_restore + f_exc + f_rad_ + f_pto + (f_moor .* n_lines) + f_damp_surge)

    @inbounds for (i, j) in enumerate(velrange(buoy_sim))
        du[j] = du1a8[i]
    end

    @inbounds for (i, j) in enumerate(posrange(buoy_sim))
        du[j] = v_vector[i]
    end
    return nothing
end

@inline function calc_mooring_force!(du, u, t, buoy_sim::FloatingBodySim{NMODES,NBODY}) where {NMODES,NBODY}
    moor_integ = buoy_sim.mooring
    r_local = buoy_sim.rloc_fixed  # assumes buoy is also at 0,0 on global
    new_time = t # current time 

    # Current main buoy position and velocity (local to CoM)
    v_buoy = get_main_buoy_vel(u, buoy_sim)
    r_buoy = get_main_buoy_pos(u, buoy_sim)

    # Precompute rotation parameters
    θ = r_buoy[end] # PARAMETRIZAR DEPOIS
    θ_dot = v_buoy[end] # PARAMETRIZAR DEPOIS
    R = rotation(θ) # PARAMETRIZAR DEPOIS

    if t > moor_integ.t
        # Update mooring solution: body CM -> fairlead position
        for i in eachindex(moor_integ.p.semis)
            r_global = rotate_point(R, r_local[i]) # R(θ) rᴮᶠ
            # println(r_global)
            ω = omega(θ_dot)

            r_buoy_ = SVector(r_buoy[1], 0.0, r_buoy[2])
            v_buoy_ = SVector(v_buoy[1], 0.0, v_buoy[2])

            r_f = point_position(r_buoy_, r_global) # assumes buoy is also at 0,0 on global at start...
            v_f = point_velocity(v_buoy_, ω, r_global) # ω = θ_dot

            input_bc = SVector(r_f..., v_f...) # INVERTER ORDEM # DEPOIS PARAMETRIZAR...
            set_boundary_conditions!(moor_integ, 1, i, new_time, input_bc)
        end
    end
    if t > moor_integ.t
        step_mooring_solver!(moor_integ)
    end

    # Step mooring solver to current time
    fmoor = get_mooring_force(moor_integ)
    fmoor_generalised = sum(generalised_force(rotate_point(R, r_local[i]), SVector(fmoor[i]...)) for i in eachindex(r_local))
    return -SVector(fmoor_generalised[1], fmoor_generalised[3], fmoor_generalised[5], 0.0, 0.0, 0.0, 0.0, 0.0)
end

@inline function calc_excitation_force!(du, u, t, buoy_sim::FloatingBodySim{NMODES}) where {NMODES}
    exc = buoy_sim.excitation

    forces = MVector{NMODES,eltype(u)}(undef)
    fill!(forces, zero(eltype(u)))

    @inbounds for k in eachindex(exc.omega)
        ωt = exc.omega[k] * t
        @simd for j = 1:NMODES
            forces[j] += exc.mag[k, j] *
                cos(ωt + exc.phase[k, j])
        end
    end

    ramp = buoy_sim.wave_spectra.ramp
    if !iszero(ramp)
        scale = min(t/ramp, one(eltype(u)))
        @inbounds @simd for j=1:NMODES
            forces[j] *= scale
        end
    end
    return SVector(forces)
end

@inline function calc_radiation_force!(du, u, t, buoy_sim::FloatingBodySim{NMODES}) where {NMODES}

    f_radiation = MVector{NMODES,eltype(u)}(undef)
    fill!(f_radiation, zero(eltype(u)))

    @inbounds for rad in buoy_sim.radiation
        state = @views u[rad.i0:rad.i1]
        dst = @views du[rad.i0:rad.i1]

        mul!(dst, rad.A, state) # ẋ = A*x
        @. dst += rad.B * u[rad.dof_j] # ẋ += B*v

        f_radiation[rad.dof_i] += dot(rad.C, state) # F = C*x
    end
    return -SVector(f_radiation)

end

@inline function calc_pto_force!(du, u, t, buoy_sim::FloatingBodySim{NMODES}) where {NMODES}
    @unpack modes, OWCs_tup, PTO_ranges = buoy_sim

    v_vector = get_buoy_vel(u, buoy_sim)
    z_vector = get_buoy_pos(u, buoy_sim)

    v1, v3, v5, v9, v15, v21, v27, v33 = v_vector # states[1:8]
    z1, z3, z5, z9, z15, z21, z27, z33 = z_vector # states[9:16]

    v_pistons = SVector{length(OWCs_tup)}(v9, v15, v21, v27, v33)
    z_pistons = SVector{length(OWCs_tup)}(z9, z15, z21, z27, z33)

    v_piston_rel = v3 .- v_pistons
    z_piston_rel = z3 .- z_pistons

    p_x_area = MVector{length(OWCs_tup),eltype(u)}(zeros(length(OWCs_tup)))
    @inbounds for i in eachindex(OWCs_tup)
        owc = OWCs_tup[i]

        du_owc = @views du[PTO_ranges[i]]
        u_owc = @views u[PTO_ranges[i]]

        rhs_pto!(du_owc, u_owc, t, z_piston_rel[i], v_piston_rel[i], owc)

        p_rel_i = extract_p_rel(u_owc, owc.turb)
        p_x_area[i] = - p_rel_i * owc.chamber_area
    end
    return SVector{NMODES,eltype(u)}(0.0, -sum(p_x_area), 0.0, p_x_area...)
end

#
end # muladd