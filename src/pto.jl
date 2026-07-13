@muladd begin # Enable FMA
#! format: noindent

# Supertypes
abstract type AbstractPTO end
abstract type AbstractTurbine end

mutable struct OWC_piston{RealT<:Real,TURB<:AbstractTurbine} <: AbstractPTO
    const idx::Int # indice do modo = IDX

    # Turbine
    turb::TURB
    const n_turbines::Int

    # Geometry
    const chamber_height::RealT
    chamber_area::RealT
end
# function OWC_piston{RealT}(idx::Int; chamber_height::Real, n_turbines::Int, turbine::AbstractTurbine) where RealT
#     println("`chamber_area` not provided")
#     return OWC_piston{RealT,typeof(turbine)}(idx, turbine, n_turbines, convert(RealT, chamber_height), convert(RealT, 0.0)) # area changed later
# end
function OWC_piston{RealT}(idx::Int; chamber_height::Real, chamber_area::Real, n_turbines::Int, turbine::AbstractTurbine) where RealT
    return OWC_piston{RealT,typeof(turbine)}(idx, turbine, n_turbines, convert(RealT, chamber_height), convert(RealT, chamber_area))
end
@inline ndofs(owc::OWC_piston) = 1

# Types of turbine
struct NoTurbine <: AbstractTurbine end

struct WellsTurbine{RealT<:Real} <: AbstractTurbine
    diam::RealT
    Psi_max::RealT
    Psi_bep::RealT
    I_turb::RealT

    pfa::RealT
    pfb::RealT

    P_rated::RealT
    function WellsTurbine(diam::RealT, pfb::Real, Psi_max::Real, Psi_bep::Real, I_turb_ref::Real, P_rated::Real, ρ_ref::Real) where RealT<:Real
        Pi_bep = Pi_Psi_bir(Psi_bep)
        pfa = ρ_ref * diam^5 * Pi_bep

        d_turb_ref = 0.5 # O QUE E ISTO?

        I_turb = I_turb_ref * (diam / d_turb_ref)^5
        new{RealT}(diam, convert(RealT, Psi_max), convert(RealT, Psi_bep), convert(RealT, I_turb), convert(RealT, pfa), convert(RealT, pfb), convert(RealT, P_rated))
    end
end

function WellsTurbine{RealT}(diam::Real; pfb, Psi_max, Psi_bep, I_turb_ref, P_rated, ρ_ref=ρ_AIR) where RealT
    return WellsTurbine(convert(RealT, diam), pfb, Psi_max, Psi_bep, I_turb_ref, P_rated, ρ_ref)
end
# WellsTurbine(diam::Real=1.2) = WellsTurbine(diam, 0.31, 0.0634) # FIX ::Real# diam= 1.2
@inline turbine_name(::WellsTurbine) = "Wells Turbine"

struct BiradialTurbine{RealT<:Real} <: AbstractTurbine
    diam::RealT
    Psi_max::RealT
    Psi_bep::RealT
    I_turb::RealT

    pfa::RealT
    pfb::RealT

    P_rated::RealT
    function BiradialTurbine(diam::RealT, pfb::Real, Psi_max::Real, Psi_bep::Real, I_turb_ref::Real, P_rated::Real, ρ_ref::Real) where RealT<:Real
        Pi_bep = Pi_Psi_bir(Psi_bep)
        pfa = ρ_ref * diam^5 * Pi_bep

        d_turb_ref = 0.5 

        I_turb = I_turb_ref * (diam / d_turb_ref)^5
        new{RealT}(diam, convert(RealT, Psi_max), convert(RealT, Psi_bep), convert(RealT, I_turb), convert(RealT, pfa), convert(RealT, pfb), convert(RealT, P_rated))
    end
end
function BiradialTurbine{RealT}(diam::Real; pfb, Psi_max, Psi_bep, I_turb_ref, P_rated, ρ_ref=ρ_AIR) where RealT
    return BiradialTurbine(convert(RealT, diam), pfb, Psi_max, Psi_bep, I_turb_ref, P_rated, ρ_ref)
end
@inline turbine_name(::BiradialTurbine) = "Biradial Turbine"

@inline ndofs(turb::Union{WellsTurbine,BiradialTurbine}) = 3
@inline extract_p_rel(u, turb::BiradialTurbine) = u[1]
@inline extract_Ω(u, turb::BiradialTurbine) = u[2]
@inline extract_E_turb(u, turb::BiradialTurbine) = u[3]

@inline function rhs_pto!(du_pto, u_pto, t, z_piston_rel, v_piston_rel, PTO::OWC_piston)
    @unpack turb, chamber_area, chamber_height, n_turbines = PTO

    # p_rel, Ω, E_turb = u_pto # State variables: p, Ω, E (3dof)
    p_rel = extract_p_rel(u_pto, turb)
    Ω = extract_Ω(u_pto, turb)
    E_turb = extract_E_turb(u_pto, turb)

    p_abs = p_rel + P_ATM
    ρ_chamber = ρ_AIR * (p_abs / P_ATM)^(1.0/γ_GAS)
    ρ_inlet = max(ρ_chamber, ρ_AIR)

    V_chamber = chamber_area * (chamber_height + z_piston_rel)
    dot_V_chamber = chamber_area * v_piston_rel
    mass_chamber = ρ_chamber * V_chamber

    Ψ = p_rel / (ρ_inlet * (Ω * turb.diam)^2)

    Φ = Phi_Psi(turb, Ψ)
    Π = Pi_Psi(turb, Ψ)

    pfa = turb.pfa
    pfb = turb.pfb

    mass_flow_single = ρ_inlet * Ω * turb.diam^3 * Φ
    mass_flow_total = mass_flow_single * n_turbines

    P_turb_single = ρ_inlet * Ω^3 * turb.diam^5 * Π
    P_gen_single = pfa * Ω^pfb

    varA = dot_V_chamber / V_chamber
    varB = mass_flow_total / mass_chamber
    denom2_1 = 1.0 / (turb.I_turb * Ω)

    du_pto[1] = -γ_GAS * p_abs * (varA + varB)
    du_pto[2] = (P_turb_single - P_gen_single) * denom2_1
    du_pto[3] = P_turb_single * n_turbines

    return nothing
end

# Functions
@inline function Phi_Psi(t::WellsTurbine, Psi)
    Psi_abs = min(abs(Psi), t.Psi_max)
    phi = if Psi_abs <= 5.247313e-03
        evalpoly(Psi_abs, (0.0, 2.15164, -133.811))
    elseif Psi_abs <= 0.09612150327465495
        evalpoly(Psi_abs, (0.00368439, 0.747339))
    else
        evalpoly(Psi_abs, (-0.0659106, 1.74814, -2.87934))
    end
    return sign(Psi) * phi
end

@inline function Pi_Psi(t::WellsTurbine, Psi)
    Psi_abs = min(abs(Psi), t.Psi_max)
    pi = if Psi_abs <= 0.11936412532777915
        evalpoly(Psi_abs, (-0.0000501325, 0.00786153, -0.332857, 18.3816, -31.2077, -2028.33, 10855.7))
    elseif Psi_abs <= 0.21291302883887459
        evalpoly(Psi_abs, (0.100773, -1.94323, 13.8623, -42.6193, 46.6723))
    else
        0.0
    end
    return pi
end

@inline function eta_Psi(t::WellsTurbine, Psi)
    Psi_abs = abs(Psi)
    num = Pi_Psi(t, Psi_abs)
    den = max(1e-4, Phi_Psi(t, Psi_abs) * Psi_abs)
    return num / den
end

@inline function Phi_Psi(t::BiradialTurbine, Psi)
    Psi_abs = abs(Psi)
    phi = (0.0 + 539581725.354 * Psi_abs + 638407403.112 * Psi_abs^2) /
          (844010065.844 + 4360645489.11 * Psi_abs + 468212474.631 * Psi_abs^2 + Psi_abs^3)
    return sign(Psi) * phi
end

@inline function Pi_Psi(t::BiradialTurbine, Psi)
    Psi_abs = abs(Psi)
    return (-8311615.56473 + 69022248.16 * Psi_abs + 136450487.222 * Psi_abs^2) /
           (1067882948.21 + 512641214.194 * Psi_abs + 76229.0830089 * Psi_abs^2 + Psi_abs^3)
end
@inline function Pi_Psi_bir(Psi)
    Psi_abs = abs(Psi)
    return (-8311615.56473 + 69022248.16 * Psi_abs + 136450487.222 * Psi_abs^2) /
           (1067882948.21 + 512641214.194 * Psi_abs + 76229.0830089 * Psi_abs^2 + Psi_abs^3)
end

@inline function eta_Psi(t::BiradialTurbine, Psi)
    Psi_abs = abs(Psi)
    return (-0.0698294618409 + 0.408116075959 * Psi_abs + 4.38745266305 * Psi_abs^2) /
           (0.073510236828 - 0.259269785174 * Psi_abs + 7.08676177897 * Psi_abs^2 + Psi_abs^3)
end


end # muladd