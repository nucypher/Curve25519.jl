#=
Ported from joined backend/serial/u64/constants.rs and constants.rs

Precomputed constants.
=#


# Edwards `2*d` value, equal to `2*(-121665/121666) mod (2^255-19)`.
const EDWARDS_D2 = convert(
    InternalScalar,
    16295367250680780974490674513165176452449235426866156013048779062215315747161)


# The Ed25519 basepoint, as an `EdwardsPoint`.
const ED25519_BASEPOINT = EdwardsPoint(convert.(
    InternalScalar, (
        15112221349535400772501151409588531511454012693041857206046113283949847762202,
        46316835694926478169428394003475163141307993866256225615783033603165251855960,
        1,
        46827403850823179245072216630277197565144205554125654976674165829533817101731,
        ))...)


# The Ristretto basepoint, as a `RistrettoPoint`.
const RISTRETTO_BASEPOINT = RistrettoPoint(ED25519_BASEPOINT)


#=
BASEPOINT_ORDER` is the order of the Ristretto group and of the Ed25519 basepoint, i.e.,
2^252 + 27742317777372353535851937790883648493.
=#
const BASEPOINT_ORDER = (one(BigInt) << 252) + 27742317777372353535851937790883648493
