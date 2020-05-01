module Curve25519

using DarkIntegers
using DarkCurves
using ConstantTime

const CT = ConstantTime

include("field.jl")

const InternalScalar = FieldElement51

include("edwards.jl")
include("internal_points.jl")
include("ristretto.jl")
include("window.jl")
include("scalar.jl")
include("constants.jl")


const RistrettoPointVT = RistrettoPoint{InternalScalar}
const RistrettoPointCT = RistrettoPoint{CT.Value{InternalScalar}}


const _ScalarType = MLUInt{2, UInt128}
const _SCALAR_MODULUS = convert(_ScalarType, BASEPOINT_ORDER)

const RistrettoScalarVT = ModUInt{_ScalarType, _SCALAR_MODULUS}

# FIXME: this assumes that MgModUInt is constant time, which is not true.
# For now we will have coarse control over constant-timeness, should be fixed
# when DarkIntegers support constant time operations.
const RistrettoScalarCT = CT.Value{RistrettoScalarVT}

Base.:>>(x::RistrettoScalarCT, s) = CT.wrap(CT.unwrap(x) >> s)
Base.rem(x::RistrettoScalarCT, tp::Type) = CT.wrap(CT.unwrap(x) % tp)



base_point(::Type{RistrettoPointVT}) = RISTRETTO_BASEPOINT
base_point(::Type{RistrettoPointCT}) = CT.wrap(RISTRETTO_BASEPOINT)


curve_order(::Type{RistrettoScalarVT}) = convert(RistrettoScalarVT, BASEPOINT_ORDER)
curve_order(::Type{RistrettoScalarCT}) = CT.wrap(curve_order(RistrettoScalarVT))

end
