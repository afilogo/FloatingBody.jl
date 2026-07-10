mutable struct NRRandom
    u::UInt64
    v::UInt64
    w::UInt64

    function NRRandom(seed::Int)
        rng = new(0, 0, 0)
        seed_rng!(rng, UInt64(seed))
        return rng
    end
end
function seed_rng!(rng::NRRandom, j::UInt64)
    rng.v = 4101842887655102017
    rng.w = 1
    rng.u = j ⊻ rng.v
    _int64!(rng)
    rng.v = rng.u
    _int64!(rng)
    rng.w = rng.v
    _int64!(rng)
    nothing
end
function _int64!(rng::NRRandom)
    # u update
    rng.u = rng.u * 2862933555777941757 + 7046029254386353087

    # v update
    rng.v = rng.v ⊻ (rng.v >> 17)
    rng.v = rng.v ⊻ (rng.v << 31)
    rng.v = rng.v ⊻ (rng.v >> 8)

    # w update
    rng.w = UInt64(4294957665) * (rng.w & 0xffffffff) + (rng.w >> 32)

    # x scrambling
    x = rng.u ⊻ (rng.u << 21)
    x = x ⊻ (x >> 35)
    x = x ⊻ (x << 4)
    return (x + rng.v) ⊻ rng.w
end

function rand_float(rng::NRRandom)
    return 5.42101086242752217e-20 * Float64(_int64!(rng))
end