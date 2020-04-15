struct FieldElement51
    limbs :: NTuple{5, UInt64}
end


Base.getindex(x::FieldElement51, i) = x.limbs[i]


Base.zero(::Type{FieldElement51}) = FieldElement51((0, 0, 0, 0, 0))


Base.one(::Type{FieldElement51}) = FieldElement51((1, 0, 0, 0, 0))


Base.:+(x::FieldElement51, y::FieldElement51) =
    FieldElement51((x[1] + y[1], x[2] + y[2], x[3] + y[3], x[4] + y[4], x[5] + y[5]))


function Base.:-(x::FieldElement51, y::FieldElement51)
    #=
    // To avoid underflow, first add a multiple of p.
    // Choose 16*p = p << 4 to be larger than 54-bit _rhs.
    //
    // If we could statically track the bitlengths of the limbs
    // of every FieldElement51, we could choose a multiple of p
    // just bigger than _rhs and avoid having to do a reduction.
    //
    // Since we don't yet have type-level integers to do this, we
    // have to add an explicit reduction call here.
    =#
    reduce(FieldElement51((
        (x[1] + 0x007ffffffffffed0) - y[1],
        (x[2] + 0x007ffffffffffff0) - y[2],
        (x[3] + 0x007ffffffffffff0) - y[3],
        (x[4] + 0x007ffffffffffff0) - y[4],
        (x[5] + 0x007ffffffffffff0) - y[5],
        )))
end


function Base.:-(x::FieldElement51)
    #// See commentary in the Sub impl
    reduce(FieldElement51((
        0x007ffffffffffed0 - x[1],
        0x007ffffffffffff0 - x[2],
        0x007ffffffffffff0 - x[3],
        0x007ffffffffffff0 - x[4],
        0x007ffffffffffff0 - x[5],
        )))
end


#/// Given 64-bit input limbs, reduce to enforce the bound 2^(51 + epsilon).
function reduce(x::FieldElement51)
    LOW_51_BIT_MASK = (one(UInt64) << 51) - one(UInt64)

    #=
    // Since the input limbs are bounded by 2^64, the biggest
    // carry-out is bounded by 2^13.
    //
    // The biggest carry-in is c4 * 19, resulting in
    //
    // 2^51 + 19*2^13 < 2^51.0000000001
    //
    // Because we don't need to canonicalize, only to reduce the
    // limb sizes, it's OK to do a "weak reduction", where we
    // compute the carry-outs in parallel.
    =#

    c1 = x[1] >> 51
    c2 = x[2] >> 51
    c3 = x[3] >> 51
    c4 = x[4] >> 51
    c5 = x[5] >> 51

    x1 = x[1] & LOW_51_BIT_MASK
    x2 = x[2] & LOW_51_BIT_MASK
    x3 = x[3] & LOW_51_BIT_MASK
    x4 = x[4] & LOW_51_BIT_MASK
    x5 = x[5] & LOW_51_BIT_MASK

    x1 += c5 * 19
    x2 += c1
    x3 += c2
    x4 += c3
    x5 += c4

    FieldElement51((x1, x2, x3, x4, x5))
end


function Base.:*(x::FieldElement51, y::FieldElement51)
    #/// Helper function to multiply two 64-bit integers with 128
    #/// bits of output.

    # Alias self, _rhs for more readable formulas
    a = x
    b = y

    #=
    // Precondition: assume input limbs a[i], b[i] are bounded as
    //
    // a[i], b[i] < 2^(51 + b)
    //
    // where b is a real parameter measuring the "bit excess" of the limbs.

    // 64-bit precomputations to avoid 128-bit multiplications.
    //
    // This fits into a u64 whenever 51 + b + lg(19) < 64.
    //
    // Since 51 + b + lg(19) < 51 + 4.25 + b
    //                       = 55.25 + b,
    // this fits if b < 8.75.
    =#
    b1_19 = b[2] * 19
    b2_19 = b[3] * 19
    b3_19 = b[4] * 19
    b4_19 = b[5] * 19

    # Multiply to get 128-bit coefficients of output
    c0 = widemul(a[1],b[1]) + widemul(a[5],b1_19) + widemul(a[4],b2_19) + widemul(a[3],b3_19) + widemul(a[2],b4_19)
    c1 = widemul(a[2],b[1]) + widemul(a[1],b[2])  + widemul(a[5],b2_19) + widemul(a[4],b3_19) + widemul(a[3],b4_19)
    c2 = widemul(a[3],b[1]) + widemul(a[2],b[2])  + widemul(a[1],b[3])  + widemul(a[5],b3_19) + widemul(a[4],b4_19)
    c3 = widemul(a[4],b[1]) + widemul(a[3],b[2])  + widemul(a[2],b[3])  + widemul(a[1],b[4])  + widemul(a[5],b4_19)
    c4 = widemul(a[5],b[1]) + widemul(a[4],b[2])  + widemul(a[3],b[3])  + widemul(a[2],b[4])  + widemul(a[1],b[5])

    #=
    // How big are the c[i]? We have
    //
    //    c[i] < 2^(102 + 2*b) * (1+i + (4-i)*19)
    //         < 2^(102 + lg(1 + 4*19) + 2*b)
    //         < 2^(108.27 + 2*b)
    //
    // The carry (c[i] >> 51) fits into a u64 when
    //    108.27 + 2*b - 51 < 64
    //    2*b < 6.73
    //    b < 3.365.
    //
    // So we require b < 3 to ensure this fits.
    =#
    #=
    debug_assert!(a[0] < (1 << 54)); debug_assert!(b[0] < (1 << 54));
    debug_assert!(a[1] < (1 << 54)); debug_assert!(b[1] < (1 << 54));
    debug_assert!(a[2] < (1 << 54)); debug_assert!(b[2] < (1 << 54));
    debug_assert!(a[3] < (1 << 54)); debug_assert!(b[3] < (1 << 54));
    debug_assert!(a[4] < (1 << 54)); debug_assert!(b[4] < (1 << 54));
    =#

    #=
    // Casting to u64 and back tells the compiler that the carry is
    // bounded by 2^64, so that the addition is a u128 + u64 rather
    // than u128 + u128.
    =#

    LOW_51_BIT_MASK = (one(UInt64) << 51) - one(UInt64)

    c1 += (c0 >> 51) % UInt64
    out0 = (c0 % UInt64) & LOW_51_BIT_MASK

    c2 += (c1 >> 51) % UInt64
    out1 = (c1 % UInt64) & LOW_51_BIT_MASK

    c3 += (c2 >> 51) % UInt64
    out2 = (c2 % UInt64) & LOW_51_BIT_MASK

    c4 += (c3 >> 51) % UInt64
    out3 = (c3 % UInt64) & LOW_51_BIT_MASK

    carry = (c4 >> 51) % UInt64
    out4 = (c4 % UInt64) & LOW_51_BIT_MASK

    #=
    // To see that this does not overflow, we need out[0] + carry * 19 < 2^64.
    //
    // c4 < a0*b4 + a1*b3 + a2*b2 + a3*b1 + a4*b0 + (carry from c3)
    //    < 5*(2^(51 + b) * 2^(51 + b)) + (carry from c3)
    //    < 2^(102 + 2*b + lg(5)) + 2^64.
    //
    // When b < 3 we get
    //
    // c4 < 2^110.33  so that carry < 2^59.33
    //
    // so that
    //
    // out[0] + carry * 19 < 2^51 + 19 * 2^59.33 < 2^63.58
    //
    // and there is no overflow.
    =#
    out0 += carry * 19

    #// Now out[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
    out1 += out0 >> 51
    out0 &= LOW_51_BIT_MASK

    #// Now out[i] < 2^(51 + epsilon) for all i.
    FieldElement51((out0, out1, out2, out3, out4))
end

