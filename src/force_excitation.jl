@muladd begin # Enable FMA
#! format: noindent

struct ExcitationData{N,RealT<:Real}
    mag::Matrix{RealT}      # Nfreq × Nmode
    phase::Matrix{RealT}    # Nfreq × Nmode
    omega::Vector{RealT}    # Nfreq
end
function ExcitationData(excData::Dict{String,Vector{T}}, omega::Vector{T}, modes::NTuple{N,Int}) where {N,T}

    nfreq = length(omega)

    mag = Matrix{T}(undef, nfreq, N)
    phase = Matrix{T}(undef, nfreq, N)

    @inbounds for (j, m) in enumerate(modes)
        mag[:, j] .= excData["$(m)_MagFe"]
        phase[:, j] .= excData["$(m)_PhaseFe"]
    end

    return ExcitationData{N,T}(mag, phase, omega)
end
@inline function excitation_force(excitation::ExcitationData, t, j)
    s = zero(eltype(excitation.omega))
    for k in eachindex(excitation.omega) # @inbounds @simd 
        s += excitation.mag[k, j] *
             cos(excitation.omega[k]*t + excitation.phase[k, j])
    end
    return s
end
function readExcitation(h5f, wave_spectra, modes, ULength)

    period = read(h5f["data/period"])
    Tbl_omega = 2.0 * π ./ period

    # Ensure Tbl_omega is sorted for interpolation
    sort_idx = sortperm(Tbl_omega)
    Tbl_omega = Tbl_omega[sort_idx]

    omega_min = minimum(Tbl_omega)
    omega_max = maximum(Tbl_omega)

    Hs = wave_spectra.Hs
    Te = wave_spectra.Te
    n_waves = wave_spectra.n_waves

    A_w, omega, RndPhase = PM_spectrum(Hs, Te, n_waves, omega_min, omega_max)

    excData = Dict{String,Vector{Float64}}()

    sf = get_scaling_factor_Xi(3, ρ_WATER, GRAV, ULength)

    for mode in modes

        Tbl_PhaseFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/phase"])[sort_idx] .* (π / 180.0)

        Tbl_MagFe = read(h5f["data/wave_heading_0.0/excitation/mode_$mode/magnitude"])[sort_idx]
        Tbl_MagFe .*= sf

        Phi = interp_vector(omega, Tbl_omega, Tbl_PhaseFe)
        PhaseFe = Phi .+ RndPhase

        Gamma = interp_vector(omega, Tbl_omega, Tbl_MagFe)
        MagFe = A_w .* Gamma

        excData["$(mode)_MagFe"] = MagFe
        excData["$(mode)_PhaseFe"] = PhaseFe

    end
    # Wave power per unit crest length (wave energy flux) for a monochromatic (single frequency), linear sinusoidal wave traveling in deep water [W/m]
    # Total power of a wave spectrum by summing up the energy of individual wave components
    unit_power_spec = sum(0.25 * ρ_WATER * (GRAV * w_a)^2 / w
                          for (w_a, w) in zip(A_w, omega)) # $P \approx 0.49 H_s^2 T_e$
    wave_power = Hs * unit_power_spec # total wave power 

    return excData, omega, A_w, wave_power
end

## Interpolation helper
function interpolate_1d(x_new, x_arr, y_arr)
    # Simple linear interpolation replacement for np.interp
    idx = searchsortedfirst(x_arr, x_new)
    if idx == 1
        ;
        return y_arr[1];
    end
    if idx > length(x_arr)
        ;
        return y_arr[end];
    end
    x0, x1 = x_arr[idx-1], x_arr[idx]
    y0, y1 = y_arr[idx-1], y_arr[idx]
    return y0 + (x_new - x0) * (y1 - y0) / (x1 - x0)
end

function interp_vector(x_new_vec, x_arr, y_arr)
    return [interpolate_1d(x, x_arr, y_arr) for x in x_new_vec]
end


#
end # muladd