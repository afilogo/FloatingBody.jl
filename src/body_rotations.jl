"""
body_rotations.jl

Only works with 3d vectors
"""

using Rotations
import LinearAlgebra: ×

@muladd begin # Enable FMA
#! format: noindent

# Cross product
×(a::SVector{2}, b::SVector{2}) = a[1]*b[2] - a[2]*b[1] # unroll for 2d

# omega vector
@inline omega(θ_dot) = SVector(zero(eltype(θ_dot)), θ_dot, zero(eltype(θ_dot)))
@inline omega(φ_dot, θ_dot, ψ_dot) = SVector(φ_dot, θ_dot, ψ_dot)

# Constructors
"Rotation matrix R(θ) for 2D (surge, heave, pitch)."
@inline rotation(θ::Real) = RotY(θ)                           # 2D # Angle2d(θ) 

"Rotation matrix for 3D Euler angles (φ, θ, ψ). Same as RotX(φ)*RotY(θ)*RotZ(ψ)"
@inline rotation(φ::Real, θ::Real, ψ::Real) = RotXYZ(φ, θ, ψ)    # 3D, swap for RotZYX etc.

"Passthrough for a pre-built rotation (e.g. QuatRotation)."
@inline rotation(q::QuatRotation) = q                            # not finished...

# Kinematics
"Body-to-global rotation of a local vector: r_g = R(θ) rᴮ."
@inline rotate_point(R, r_local) = SVector(R * r_local)

"Fairlead position in the absolute frame: rᶠ = r_o + R(θ) rᴮᶠ."
@inline point_position(r, r_global) = SVector(r + r_global)

"Fairlead velocity in the absolute frame: vᶠ = v_o + ω × r_global"
@inline point_velocity(v, ω, r_global) = SVector(v + ω × r_global)


# Generalised mooring force

"""
Generalised force Qmoor = Jᵀ Fmoor
- 2D: [Fx, Fz, My],  My = r_g × Fmoor 
- 3D: [Fx, Fy, Fz, Mx, My, Mz]

r_g = R(θ) rᴮᶠ : fairlead offset from buoy origin, expressed in the global frame.
f             : mooring force in the global frame.
"""
# @inline generalised_force(r_global::SVector{2}, f::SVector{2}) = SVector(f[1], f[2], r_global[1]*f[2] - r_global[2]*f[1])
@inline generalised_force(r_global::SVector{2}, f::SVector{3}) = SVector(f[1], f[2], r_global[1]*f[3] - r_global[2]*f[1])
@inline generalised_force(r_global::SVector{3}, f::SVector{3}) = SVector(f[1], f[2], f[3], (r_global × f)...)
#
# USAGE:
# R        = rotation(θ)
# r_global = rotate_point(R, r_local) # R(θ) rᴮᶠ — used everywhere below

# r_f = point_position(r, r_global)
# v_f = point_velocity(v, ω, r_global)
# Qmoor = generalised_force(r_global, f)
end # muladd 