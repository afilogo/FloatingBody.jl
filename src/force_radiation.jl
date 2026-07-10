@muladd begin # Enable FMA
#! format: noindent

# CALCULATE FORCE WITH PRONY'S METHOD.

struct RadiationSS{T}
    dof_i::Int
    dof_j::Int

    A::Matrix{T}
    B::Vector{T}
    C::Vector{T}

    i0::Int
    i1::Int
end
function readRadiation(filename::String, modes::AbstractVector{<:Integer}, next_state::Int)
    radiation = RadiationSS{Float64}[]

    # next_state = 31

    h5open(filename, "r") do h5f
        for mode in sort!(collect(keys(h5f)))

            dof_i = Int(read(h5f["$mode/dof_i"]))
            dof_j = Int(read(h5f["$mode/dof_j"]))

            idx_i = findfirst(==(dof_i), modes)
            idx_j = findfirst(==(dof_j), modes)

            idx_i === nothing && error("Mode $dof_i not found in modes.")
            idx_j === nothing && error("Mode $dof_j not found in modes.")

            A = Matrix(read(h5f["$mode/A"])')
            B = vec(Matrix(read(h5f["$mode/B"])'))
            C = vec(Matrix(read(h5f["$mode/C"])'))

            n = size(A, 1)

            i0 = next_state
            i1 = next_state + n - 1

            push!(radiation,
                RadiationSS(
                    idx_i,
                    idx_j,
                    A,
                    B,
                    C,
                    i0,
                    i1,
                )
            )

            next_state = i1 + 1
        end
    end

    return radiation, next_state

end
#
end # muladd