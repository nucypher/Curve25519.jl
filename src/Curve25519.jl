module Curve25519

using DarkIntegers


SType = MLUInt{2, UInt128}
InternalMod = convert(SType, (big(1) << 255) - 19)
ModSType = MgModUInt{SType, InternalMod}
RistrettoScalar = ModSType


function from_FieldElement51(x::Array{T, 1}) where T
    @assert length(x) == 5
    xx = big.(x)
    convert(ModSType, xx[1] + (xx[2] << 51) + (xx[3] << 102) + (xx[4] << 153) + (xx[5] << 204))
end


function from_bytes(x::Array{T, 1}) where T
    @assert length(x) == 32
    res = zero(BigInt)
    for i in 1:32
        res += big(x[i]) << ((i - 1) * 8)
    end
    convert(ModSType, res)
end


include("field.jl")
include("montgomery.jl")
include("edwards.jl")
include("ristretto.jl")
include("window.jl")
include("scalar.jl")

include("backend/serial/curve_models.jl")
include("backend/serial/scalar_mul/variable_base.jl")

include("backend/serial/constants.jl")
include("constants.jl")

base_point() = RISTRETTO_BASEPOINT_POINT
curve_order() = BASEPOINT_ORDER

end
