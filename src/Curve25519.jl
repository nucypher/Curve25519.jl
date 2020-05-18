module Curve25519

using DarkIntegers
using DarkCurves
using ConstantTime
using Random

const CT = ConstantTime

include("field.jl")

const InternalScalar = FieldElement51

include("edwards.jl")
include("internal_points.jl")
include("ristretto.jl")
include("window.jl")
include("scalar.jl")
include("constants.jl")

export RistrettoCurveVT

const _ScalarType = MLUInt{4, UInt64}
const _SCALAR_MODULUS = convert(_ScalarType, BASEPOINT_ORDER)

DarkCurves.curve_point_type(::Type{RistrettoCurveVT}) = RistrettoPoint{InternalScalar}
DarkCurves.curve_scalar_type(::Type{RistrettoCurveVT}) = ModUInt{_ScalarType, _SCALAR_MODULUS}

Base.one(::Type{RistrettoPoint{InternalScalar}}) = RISTRETTO_BASEPOINT


const RistrettoScalarVT = curve_scalar_type(RistrettoCurveVT)
const RistrettoPointVT = curve_point_type(RistrettoCurveVT)


const _BASE_POWERS_TABLE = RistrettoBasepointTable(one(RistrettoPointVT))

Random.Sampler(RNG::Type{<:AbstractRNG}, ::Type{RistrettoPointVT}, n::Union{Val{1}, Val{Inf}}) =
    Random.SamplerType{RistrettoPointVT}()

function Base.rand(rng::AbstractRNG, ::Random.SamplerType{RistrettoPointVT})
    _BASE_POWERS_TABLE * rand(rng, RistrettoScalarVT)
end


end
