## Sea State Tables
module WaveClimateAtlantic
const HS_TABLE = (1.10, 1.18, 1.23, 1.88, 1.96, 2.07, 2.14, 3.06, 3.18, 3.29, 4.75, 4.91)
const TE_TABLE = (5.49, 6.5, 7.75, 6.33, 7.97, 9.75, 11.58, 8.03, 9.93, 11.80, 9.84, 12.03)
const TP_TABLE = (6.076, 7.193, 8.577, 7.005, 8.82, 10.790, 12.820, 8.887, 10.99, 13.060, 10.89, 13.310)
const PROB_TABLE = (7.14, 12.53, 8.29, 11.74, 20.96, 8.73, 0.60, 9.55, 10.22, 2.61, 4.79, 2.85)
end

@muladd begin # Enable FMA
#! format: noindent

struct WaveSpectra{RealT<:Real} # Irregular
    Hs::RealT
    Tp::RealT
    Te::RealT
    n_waves::Int

    ramp::RealT # ramp-up simulations
end
function WaveSpectra{RealT}(Hs::Real, Tp; n_waves::Int=400, ramp::Real=20.0) where RealT
    Te = Tp / 1.14
    return WaveSpectra{RealT}(convert(RealT, Hs), convert(RealT, Tp), convert(RealT, Te), n_waves, convert(RealT, ramp))
end

function PM_spectrum(Hs, Te, nwaves, ωmin, ωmax)
    ωspan = ωmax - ωmin
    Δω = ωspan / (nwaves - 1)

    rng = NRRandom(17)

    Δωrand = Vector{Float64}(undef, nwaves)

    sumΔω = 0.0
    @inbounds for i in eachindex(Δωrand)
        sign = rand_float(rng) <= 0.5 ? 0.2 : -0.2
        Δωrand[i] = (1 + sign * rand_float(rng)) * Δω

        sumΔω += Δωrand[i]
    end

    sumΔω -= (Δωrand[1] + Δωrand[end]) / 2
    Δωrand .*= (0.999999 * ωspan) / sumΔω

    ω = Vector{Float64}(undef, nwaves)
    amp = Vector{Float64}(undef, nwaves)
    phase = Vector{Float64}(undef, nwaves)

    ω[1] = ωmin
    @inbounds begin
        φ = 2π * rand_float(rng) # 1st frequency
        phase[1] = φ

        x = Te * ω[1]
        x4 = (x*x)*(x*x)
        S = 262.9 / (ω[1] * x4) * exp(-1054.0 / x4)
        amp[1] = Hs * sqrt(2 * Δωrand[1] * S)

        for i = 2:nwaves

            ω[i] = ω[i-1] + 0.5*(Δωrand[i] + Δωrand[i-1])

            φ = 2π * rand_float(rng)
            phase[i] = φ

            x = Te * ω[i]
            x4 = (x*x)*(x*x)
            S = 262.9 / (ω[i] * x4) * exp(-1054.0 / x4)

            amp[i] = Hs * sqrt(2 * Δωrand[i] * S)

        end
    end
    return amp, ω, phase
end

# Regular wave power
# $$\omega^2 = g k \tanh(kh)$$
@inline function dispersion_relation(ω, h)
    k0h = h * ω^2 / GRAV # deep water 
    return k0h * sqrt(1.0 + 1.0 / (1E-12 + k0h * (1.0 + 0.6522k0h + 0.4622k0h^2 + 0.0864k0h^4 + 0.0675k0h^5))) # correction
end
@inline function group_velocity(h, H, ω)
    two_kh = 2.0 * dispersion_relation(ω, h)
    c_g = ω * h / two_kh * (1.0 + two_kh / sinh(two_kh))
    return c_g
end
@inline function reg_wave_power(rho_x_g, h, H, ω)
    c_g = group_velocity(h, H, ω)
    E = 1.0 / 8.0 * rho_x_g * H^2
    return E * c_g
end


##
# filename = "results_regular" 
# config["spectrum"] = false
# config["omega"] = 0.9624993711202005
# config["Aw"] = 2.5
# config["depth"] = 40.0


# struct RegularWaves{RealT<:Real}
#     filename::String
#     modes::Vector{RealT}
#     ramp::RealT
#     # itp_1d::ITP

#     ω::RealT # omega
#     A_w::RealT
#     h_depth::RealT

#     fe_omega::RealT # ω
#     # W_amp::RealT # A_w
#     PhaseFe::Vector{RealT}
#     MagFe::Vector{RealT}
#     wave_power::RealT
# end

# struct WaveSpectra{RealT<:Real} # Irregular
#     filename::String
#     modes::Vector{RealT}
#     ramp::RealT

#     fe_omega::Vector{RealT} # ω
#     W_amp::Vector{RealT}
#     wave_power::RealT

#     PhaseFe::Matrix{RealT}
#     MagFe::Matrix{RealT} # each column is 400 waves

# end



# function IrregularWaves(wamit_post_file::String, sea_state_index::Int, n_waves::Int, wave_climate::WaveClimateAtlantic, RealT)
#     filename = "results_SS$(lpad(sea_state_index, 2, '0'))"
#     modes = [1, 3, 5, 9, 15, 21, 27, 33]

#     Hs = wave_climate.HS_TABLE[sea_state_index]
#     Tp = wave_climate.TP_TABLE[sea_state_index]
#     Te = Tp / 1.14

#     h5open(wamit_post_file, "r") do h5f # "./OCTAPLAT.h5"
#         ULength = read(h5f["data/ULength"])

#         period = read(h5f["data/period"])
#         Tbl_omega = 2.0 * π ./ period

#         # Ensure Tbl_omega is sorted for interpolation
#         sort_idx = sortperm(Tbl_omega)
#         Tbl_omega_sort = Tbl_omega[sort_idx]

#         ω_min = minimum(Tbl_omega_sort)
#         ω_max = maximum(Tbl_omega_sort)

#         W_amp, ω, RndPhase = SpectrumData(Hs, Te, n_waves, ω_min, ω_max) # size of each vector: n_waves

#         MagFe = zeros(RealT, n_waves, length(modes))
#         PhaseFe = zeros(RealT, n_waves, length(modes))

#         for (i, mode) in enumerate(modes)
#             Tbl_PhaseFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/phase"])[sort_idx] .* (π / 180.0)
#             Tbl_MagFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/magnitude"])[sort_idx]

#             sf = get_scaling_factor_Xi(mode, ρ_WATER, GRAV, ULength)
#             Tbl_MagFe = Tbl_MagFe .* sf # scale

#             Φ = [interpolate_1d(w, Tbl_omega_sort, Tbl_PhaseFe) for w in ω]
#             PhaseFe_wavei = Φ .+ RndPhase

#             Γ = [interpolate_1d(w, Tbl_omega_sort, Tbl_MagFe) for w in ω]
#             MagFe_wavei = W_amp .* Γ

#             PhaseFe[:, i] .= PhaseFe_wavei
#             MagFe[:, i] .= MagFe_wavei
#         end

#         unit_power_spec = sum((0.25 * ρ_WATER * (GRAV * w_a)^2 / w) for (w_a, w) in zip(W_amp, ω))
#         wave_power = Hs * unit_power_spec
#         fe_omega = ω

#     end
#     ramp = 20.0
#     return IrregularWaves{eltype(RealT)}(filename, modes, ramp, fe_omega, W_amp, wave_power, PhaseFe, MagFe)
# end



# function RegularWaves(wamit_post_file::String, ω::Real, A_w::Real, h_depth::Real, RealT)
#     filename = "results_regular"
#     modes = [1, 3, 5, 9, 15, 21, 27, 33]


#     ω = 0.9624993711202005
#     A_w = 2.5 # [m] 
#     h_depth = 40.0 # [m]

#     h5open(wamit_post_file, "r") do h5f # "./OCTAPLAT.h5"
#         ULength = read(h5f["data/ULength"])

#         period = read(h5f["data/period"])
#         Tbl_omega = 2.0 * π ./ period

#         # Ensure Tbl_omega is sorted for interpolation
#         sort_idx = sortperm(Tbl_omega)
#         Tbl_omega_sort = Tbl_omega[sort_idx]

#         ω_min = minimum(Tbl_omega_sort)
#         ω_max = maximum(Tbl_omega_sort)

#         MagFe = zeros(RealT, length(modes))
#         PhaseFe = zeros(RealT, length(modes))
#         for (i, mode) in enumerate(modes)
#             Tbl_PhaseFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/phase"])[sort_idx] .* (π / 180.0)
#             Tbl_MagFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/magnitude"])[sort_idx]

#             sf = get_scaling_factor_Xi(mode, ρ_WATER, GRAV, ULength)
#             Tbl_MagFe = Tbl_MagFe .* sf # scale

#             Φ = interpolate_1d(ω, Tbl_omega_sort, Tbl_PhaseFe)
#             Γ = interpolate_1d(ω, Tbl_omega_sort, Tbl_MagFe)

#             PhaseFe[i] = Φ
#             MagFe[i] = A_w * Γ
#         end


#         fe_omega = ω
#         W_amp = A_w
#         wave_power = reg_wave_power(ρ_WATER * GRAV, h_depth, A_w*2, ω)
#     end
#     ramp = 20.0
#     return RegularWaves{eltype(RealT)}(filename, modes, ramp, ω, A_w, h_depth, fe_omega, W_amp, PhaseFe, MagFe, wave_power)
# end


end # muladd