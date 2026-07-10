struct SimulationOutputs{T<:Real,NMODES,NOWCs}
    time::Vector{T}

    position::Matrix{T}
    velocity::Matrix{T}

    p_rel::Matrix{T}
    Ω::Matrix{T}
    E_turb::Matrix{T}

    Q_turb::Matrix{T}
    Ψ::Matrix{T}
    Φ::Matrix{T}
    Π::Matrix{T}
    η::Matrix{T}

    P_pneu::Matrix{T}
    P_turb::Matrix{T}
    P_gen::Matrix{T}

    F_exc::Matrix{T}
end
function SimulationOutputs{T}(nt::Int, NMODES::Int, NOWCs::Int) where T<:Real
    SimulationOutputs{T,NMODES,NOWCs}(
        zeros(T, nt), zeros(T, nt, NMODES),
        zeros(T, nt, NMODES), zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs), zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs), zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs),
        zeros(T, nt, NOWCs), zeros(T, nt, NMODES)
    )
end

function compute_output!(out, sol::SciMLBase.ODESolution, buoy_sim::FloatingBodySim{NMODES,NBODY,NOWCs,NVARS}) where {NMODES,NBODY,NOWCs,NVARS}
    @unpack OWCs_tup = buoy_sim
    out.time .= sol.t

    for t in eachindex(sol.u)
        u = sol.u[t]
        tim = sol.t[t]

        v = get_buoy_vel(u, buoy_sim)
        z = get_buoy_pos(u, buoy_sim)

        p_rel = MVector{NOWCs,eltype(u)}(zeros(NOWCs))
        Ω = similar(p_rel)
        E_turb = similar(p_rel)
        for i in eachindex(OWCs_tup) # @inbounds
            owc = OWCs_tup[i]
            u_owc = @views u[buoy_sim.PTO_ranges[i]]

            p_rel[i] = extract_p_rel(u_owc, owc.turb)
            Ω[i] = extract_Ω(u_owc, owc.turb)
            E_turb[i] = extract_E_turb(u_owc, owc.turb)
        end

        p_abs = @. p_rel+P_ATM
        ρ_chamber = @. ρ_AIR*(p_abs/P_ATM)^(1.0/γ_GAS)
        ρ_inlet = max.(ρ_chamber, ρ_AIR)

        turbs = ntuple(i -> OWCs_tup[i].turb, length(OWCs_tup))
        turbs_diam = ntuple(i -> turbs[i].diam, length(OWCs_tup))
        turbs_pfa = ntuple(i -> turbs[i].pfa, length(OWCs_tup))
        turbs_pfb = ntuple(i -> turbs[i].pfb, length(OWCs_tup))

        owcs_nturb = ntuple(i -> OWCs_tup[i].n_turbines, length(OWCs_tup))

        Ψ = @. p_rel/(ρ_inlet*(Ω^2*turbs_diam^2))

        Φ = similar(Ψ)
        Π = similar(Ψ)
        η = similar(Ψ)

        @inbounds for k in eachindex(Ψ)
            Φ[k] = Phi_Psi(turbs[k], Ψ[k])
            Π[k] = Pi_Psi(turbs[k], Ψ[k])
            η[k] = eta_Psi(turbs[k], Ψ[k])
        end

        Q_turb = @. owcs_nturb * Ω * turbs_diam^3 * Φ
        P_pneu = @. p_rel*Q_turb
        P_turb = @. owcs_nturb*ρ_inlet*Ω^3*turbs_diam^5*Π
        P_gen = @. owcs_nturb*turbs_pfa*Ω^turbs_pfb

        exc = buoy_sim.excitation

        f_exc = zeros(eltype(u), NMODES)
        @inbounds for m in 1:NMODES
            s = zero(eltype(u))
            for k in eachindex(exc.omega)
                s += exc.mag[k, m] * cos(exc.omega[k]*t + exc.phase[k, m])
            end
            f_exc[m] = s
        end

        ramp = buoy_sim.wave_spectra.ramp
        if !iszero(ramp)
            scale = min(t/ramp, one(eltype(u)))
            @inbounds @simd for j=1:NMODES
                f_exc[j] *= scale
            end
        end

        out.position[t, :] .= z
        out.velocity[t, :] .= v

        out.p_rel[t, :] .= p_rel
        out.Ω[t, :] .= Ω
        out.E_turb[t, :] .= E_turb

        out.Q_turb[t, :] .= Q_turb
        out.Ψ[t, :] .= Ψ
        out.Φ[t, :] .= Φ
        out.Π[t, :] .= Π
        out.η[t, :] .= η

        out.P_pneu[t, :] .= P_pneu
        out.P_turb[t, :] .= P_turb
        out.P_gen[t, :] .= P_gen

        out.F_exc[t, :] .= f_exc
    end
    return nothing
end