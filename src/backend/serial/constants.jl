#/// `APLUS2_OVER_FOUR` is (A+2)/4. (This is used internally within the Montgomery ladder.)
#pub(crate) const APLUS2_OVER_FOUR: FieldElement51 = FieldElement51([121666, 0, 0, 0, 0]);
APLUS2_OVER_FOUR = from_FieldElement51([121666, 0, 0, 0, 0])


#/// Edwards `d` value, equal to `-121665/121666 mod p`.
EDWARDS_D = from_FieldElement51([
    929955233495203,
    466365720129213,
    1662059464998953,
    2033849074728123,
    1442794654840575,
])
#EDWARDS_D = mod((MODULUS1 - 121665) * inv(121666, MODULUS1), MODULUS1)


# Edwards `2*d` value, equal to `2*(-121665/121666) mod p`.
EDWARDS_D2 = from_FieldElement51([
    1859910466990425,
    932731440258426,
    1072319116312658,
    1815898335770999,
    633789495995903,
])


#/// Precomputed value of one of the square roots of -1 (mod p)
SQRT_M1 = from_FieldElement51([
    1718705420411056,
    234908883556509,
    2233514472574048,
    2117202627021982,
    765476049583133,
])
#SQRT_M1 = 19681161376707505956807079304988542015446066515923890162744021073123829784752



#=
/// The Ed25519 basepoint, as an `EdwardsPoint`.
///
/// This is called `_POINT` to distinguish it from
/// `ED25519_BASEPOINT_TABLE`, which should be used for scalar
/// multiplication (it's much faster).
=#
ED25519_BASEPOINT_POINT = EdwardsPoint(
    from_FieldElement51([
        1738742601995546,
        1146398526822698,
        2070867633025821,
        562264141797630,
        587772402128613,
    ]),
    from_FieldElement51([
        1801439850948184,
        1351079888211148,
        450359962737049,
        900719925474099,
        1801439850948198,
    ]),
    from_FieldElement51([1, 0, 0, 0, 0]),
    from_FieldElement51([
        1841354044333475,
        16398895984059,
        755974180946558,
        900171276175154,
        1821297809914039,
    ]),
    )
