@muladd begin # Enable FMA
#! format: noindent

mutable struct FloatingBodySim{NMODES,NBODY,NOWCs,NVARS,MAT,OWC<:AbstractPTO,EXC<:ExcitationData,WSPEC,MOOR,R,RealT<:Real} # MOOR<:OrdinaryDiffEqCore.ODEIntegrator
    ULength::RealT # characteristic length of WAMIT simulations --> 1
    modes::NTuple{NMODES,Int}
    # PTO
    OWCs_tup::NTuple{NOWCs,OWC}
    OWC_ranges::NTuple{NOWCs,UnitRange{Int}}
    PTO_ranges::NTuple{NOWCs,UnitRange{Int}}
    chamber_areas::Vector{RealT}
    # Mass matrix of all bodies
    M::Matrix{RealT}
    Mg::Matrix{RealT}
    Mg_Inv::SMatrix{NMODES,NMODES,RealT,MAT} # inverse
    # Restoring matrix
    C::SMatrix{NMODES,NMODES,RealT,MAT}
    # Force radiation 
    radiation::Vector{RadiationSS{RealT}}
    # Wave spectra
    wave_spectra::WSPEC
    # Force excitation
    excitation::EXC
    # Mooring force
    mooring::MOOR
    rloc_fixed::R
    n_lines::RealT
    # others
    A_w::Vector{RealT}
    wave_power::RealT
    # B11::RealT # temporary mooring ELIMINAR DEPOIS
    function FloatingBodySim{RealT}(modes::AbstractVector, vec_PTOs::AbstractVector;
        body_file::String, rad_file::String, wave_spectra::WaveSpectra, mooring::Union{Nothing,ODEIntegrator}=nothing, n_lines::Real) where RealT

        modes_tup = tuple(modes...)

        N_OWCs = length(vec_PTOs)
        N_MODES = length(modes)

        body_modes = filter(<=(6), modes)
        N_BODY = length(body_modes)

        owc_modes = filter(>(6), modes)

        offset = N_MODES*2+1
        OWC_ranges = map(vec_PTOs) do owc
            r = offset:(offset+ndofs(owc)-1)
            offset += ndofs(owc)
            r
        end

        PTO_ranges = ntuple(i -> (N_MODES*2+(i-1)*ndofs(vec_PTOs[i].turb)+1):(N_MODES*2+i*ndofs(vec_PTOs[i].turb)), N_OWCs)

        # Check list of OWCs and modes
        for i in eachindex(vec_PTOs)
            owc = vec_PTOs[i]
            indice = owc.idx
            indice ∈ owc_modes || throw(ArgumentError("Erro"))
        end

        h5open(body_file, "r") do h5f
            ULength = read(h5f["data/ULength"])

            # Get the total number of bodies
            num_bodies = length(keys(h5f["data/disp_vol"]))
            num_owcs = num_bodies-1
            @assert N_OWCs == num_owcs "Number of defined OWCs does not match number of chambers found."

            println("Total number of bodies found: ", num_bodies)
            println("Total number of OWCs found: ", num_owcs)

            # Main body (body_0)
            cg_body_0 = Vector{RealT}(read(h5f["data/CG/body_0"]) .* ULength)
            cg_body_0[3] -= 10.0

            sf3 = get_scaling_factor_Aij(3, 3, ρ_WATER, ULength)
            buoy_mass = read(h5f["data/disp_vol/body_0"]) * sf3
            piston_mass = read(h5f["data/disp_vol/body_1"]) * sf3 # mass_9
            mass_15 = read(h5f["data/disp_vol/body_2"]) * sf3
            mass_21 = read(h5f["data/disp_vol/body_3"]) * sf3
            mass_27 = read(h5f["data/disp_vol/body_4"]) * sf3
            mass_33 = read(h5f["data/disp_vol/body_5"]) * sf3
            body_masses = Vector{RealT}([buoy_mass, piston_mass, mass_15, mass_21, mass_27, mass_33])
            @assert length(body_masses) == num_bodies

            A9 = read(h5f["data/restoring/body_1/C_3_3"]) * ULength^2
            A15 = read(h5f["data/restoring/body_2/C_3_3"]) * ULength^2
            A21 = read(h5f["data/restoring/body_3/C_3_3"]) * ULength^2
            A27 = read(h5f["data/restoring/body_4/C_3_3"]) * ULength^2
            A33 = read(h5f["data/restoring/body_5/C_3_3"]) * ULength^2
            chamber_areas = Vector{RealT}([A9, A15, A21, A27, A33])
            @assert length(chamber_areas) == num_owcs

            for i = 1:num_owcs
                # println(vec_PTOs[i].chamber_area)
                # println(chamber_areas[i])
                @assert vec_PTOs[i].chamber_area == chamber_areas[i]
            end
            # for i = 1:num_owcs
            #     vec_PTOs[i].chamber_area = chamber_areas[i]
            # end

            C33 = read(h5f["data/restoring/body_0/C_3_3"]) * ρ_WATER * GRAV * ULength^2
            C35 = read(h5f["data/restoring/body_0/C_3_5"]) * ρ_WATER * GRAV * ULength^4
            C53 = read(h5f["data/restoring/body_0/C_3_5"]) * ρ_WATER * GRAV * ULength^4
            C55 = read(h5f["data/restoring/body_0/C_5_5"]) * ρ_WATER * GRAV * ULength^4 - buoy_mass * GRAV * cg_body_0[3]
            C99 = read(h5f["data/restoring/body_1/C_3_3"]) * ρ_WATER * GRAV * ULength^2
            C1515 = read(h5f["data/restoring/body_2/C_3_3"]) * ρ_WATER * GRAV * ULength^2
            C2121 = read(h5f["data/restoring/body_3/C_3_3"]) * ρ_WATER * GRAV * ULength^2
            C2727 = read(h5f["data/restoring/body_4/C_3_3"]) * ρ_WATER * GRAV * ULength^2
            C3333 = read(h5f["data/restoring/body_5/C_3_3"]) * ρ_WATER * GRAV * ULength^2

            # Calculate moment of inertia
            C33_total = 3.4300
            C33_piston = 1.1549
            D_out_sq = (4.0 * C33_total * ULength^2) / π
            D_in_sq = (4.0 * C33_piston * ULength^2) / π
            Iyy = ρ_WATER * (π / 4.0) * (D_out_sq^2 - D_in_sq^2)

            # Build mass matrix M
            M = zeros(RealT, N_MODES, N_MODES)
            M[1, 1] = buoy_mass
            M[1, 3] = cg_body_0[3] * buoy_mass
            M[2, 2] = buoy_mass
            M[2, 3] = -cg_body_0[1] * buoy_mass
            M[3, 1] = cg_body_0[3] * buoy_mass
            M[3, 2] = -cg_body_0[1] * buoy_mass
            M[3, 3] = Iyy
            M[4, 4] = piston_mass
            M[5, 5] = mass_15
            M[6, 6] = mass_21
            M[7, 7] = mass_27
            M[8, 8] = mass_33

            Ainf = zeros(RealT, N_MODES, N_MODES)
            for (i, mi) in enumerate(modes)
                for (j, mj) in enumerate(modes)
                    try
                        path = "data/mode_$(mi)_$(mj)/added_mass"
                        val = read(h5f[path])[1]
                        Ainf[i, j] = val * get_scaling_factor_Aij(mi, mj, ρ_WATER, ULength)
                    catch
                        # If mode doesn't exist, keep 0.0
                    end
                end
            end
            Mg = M .+ Ainf

            C11 = 0.1 * C33
            # Mooring
            if isnothing(mooring)
                B11 = 2.0 * 0.6 * sqrt(C11 * Mg[1, 1])
                println("TD B11 = $(B11), C11 = $(C11)")
                rloc_fixed = nothing

                mooring_ = B11
            else
                semis = mooring.p.semis
                N_MOOR = length(semis) # number of moorings
                println("Mooring system found: $N_MOOR lines")

                if N_MOOR == 1
                    rloc_fixed = ntuple(i -> Vector{RealT}(semis[i].initial_condition.initial_position_top), N_MOOR)
                else
                    rloc_fixed = ntuple(i -> Vector{RealT}(semis[i].semis[1].initial_condition.initial_position_top), N_MOOR)
                end
                mooring_ = mooring
            end

            # Build restoring matrix C
            C = zeros(RealT, N_MODES, N_MODES) # PARAMETRIZAR
            # C[1, 1] = C11
            C[2, 2] = C33
            C[2, 3] = C35
            C[3, 2] = C53
            C[3, 3] = C55
            C[4, 4] = C99
            C[5, 5] = C1515
            C[6, 6] = C2121
            C[7, 7] = C2727
            C[8, 8] = C3333
            C = SMatrix{N_MODES,N_MODES}(C)

            Mg_Inv = inv(Mg)
            Mg_Inv = SMatrix{N_MODES,N_MODES}(Mg_Inv)

            next_state = PTO_ranges[end][end]
            # next_state = 31
            radiation, N_VARS = readRadiation(rad_file, modes, next_state)


            # n_outputs = 79

            excData, omega, A_w, wave_power = readExcitation(h5f, wave_spectra, modes, ULength,)

            excitation = ExcitationData(excData, omega, modes_tup)
            # excitation = ExcitationData(excData, omega, modes_tup)

            unit_power_spec = sum((0.25 * ρ_WATER * (GRAV * w_a)^2 / w) for (w_a, w) in zip(A_w, omega))
            wave_power = wave_spectra.Hs * unit_power_spec

            # Assume que so uso o mesmo tipo de owc em todas as chambers...
            return new{N_MODES,N_BODY,N_OWCs,N_VARS,N_MODES*N_MODES,typeof(vec_PTOs[1]),typeof(excitation),typeof(wave_spectra),typeof(mooring_),typeof(rloc_fixed),RealT}(
                ULength, modes_tup,
                tuple(vec_PTOs...), tuple(OWC_ranges...), tuple(PTO_ranges...), chamber_areas,
                M, Mg, Mg_Inv,
                C,
                radiation,
                wave_spectra, excitation,
                mooring_, rloc_fixed, n_lines,
                A_w, wave_power)
        end
    end
end
@inline nmodes(::FloatingBodySim{NMODES,NBODY,NOWCs,NVARS}) where {NMODES,NBODY,NOWCs,NVARS} = NMODES
@inline nowcs(::FloatingBodySim{NMODES,NBODY,NOWCs,NVARS}) where {NMODES,NBODY,NOWCs,NVARS} = NOWCs
@inline nvars(::FloatingBodySim{NMODES,NBODY,NOWCs,NVARS}) where {NMODES,NBODY,NOWCs,NVARS} = NVARS

# @inline get_buoy_vel(u, buoy_sim::FloatingBodySim{NMODES}) where NMODES = SVector{NMODES}(@views u[1:NMODES])
# @inline get_buoy_pos(u, buoy_sim::FloatingBodySim{NMODES}) where NMODES = SVector{NMODES}(@views u[(NMODES+1):(NMODES*2)])

@inline get_buoy_vel(u, buoy_sim::FloatingBodySim{NMODES}) where {NMODES} =
    SVector{NMODES}(@view u[velrange(buoy_sim)])
@inline get_buoy_pos(u, buoy_sim::FloatingBodySim{NMODES}) where {NMODES} =
    SVector{NMODES}(@view u[posrange(buoy_sim)])

@inline get_main_buoy_vel(u, buoy_sim::FloatingBodySim{NMODES,NBODY}) where {NMODES,NBODY} =
    SVector{NBODY}(@view u[velrange_mainbody(buoy_sim)])
@inline get_main_buoy_pos(u, buoy_sim::FloatingBodySim{NMODES,NBODY}) where {NMODES,NBODY} =
    SVector{NBODY}(@view u[posrange_mainbody(buoy_sim)])

@inline velrange(::FloatingBodySim{NMODES}) where {NMODES} = Base.OneTo(NMODES)
@inline posrange(::FloatingBodySim{NMODES}) where {NMODES} = (NMODES+1):(2NMODES)

@inline velrange_mainbody(::FloatingBodySim{NMODES,NBODY}) where {NMODES,NBODY} = Base.OneTo(NBODY)
@inline posrange_mainbody(::FloatingBodySim{NMODES,NBODY}) where {NMODES,NBODY} = (NMODES+1):(NMODES+NBODY)



#
end # muladd