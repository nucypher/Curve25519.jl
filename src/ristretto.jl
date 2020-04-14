struct CompressedRistretto{T}
    val :: T
end


Base.zero(::Type{CompressedRistretto{T}}) where T = CompressedRistretto{T}(zero(T))


#=
/// Attempt to decompress to an `RistrettoPoint`.
///
/// # Return
///
/// - `Some(RistrettoPoint)` if `self` was the canonical encoding of a point;
///
/// - `None` if `self` was not the canonical encoding of a point.
=#
function decompress(p::CompressedRistretto{T}) where T
    #=
    // Step 1. Check s for validity:
    // 1.a) s must be 32 bytes (we get this from the type system)
    // 1.b) s < p
    // 1.c) s is nonnegative
    //
    // Our decoding routine ignores the high bit, so the only
    // possible failure for 1.b) is if someone encodes s in 0..18
    // as s+p in 2^255-19..2^255-1.  We can check this by
    // converting back to bytes, and checking that we get the
    // original input, since our encoding routine is canonical.
    =#

    s = p.val

    # TODO: apparently checks for cases where `val` is greater than the modulus?
    # s_encoding_is_canonical = &s_bytes_check[..].ct_eq(self.as_bytes());

    s_is_negative = is_negative(s)

    if s_is_negative
        return nothing
    end

    # Step 2.  Compute (X:Y:Z:T).
    one_ = one(T)
    ss = square(s)
    u1 = one_ - ss #      //  1 + as²
    u2 = one_ + ss #      //  1 - as²    where a=-1
    u2_sqr = square(u2)  # (1 - as²)²

    # v == ad(1+as²)² - (1-as²)²            where d=-121665/121666
    v = (-EDWARDS_D) * square(u1) - u2_sqr

    (ok, I) = invsqrt(v * u2_sqr) # 1/sqrt(v*u_2²)

    Dx = I * u2 # 1/sqrt(v)
    Dy = I * (Dx * v) # 1/u2

    # x == | 2s/sqrt(v) | == + sqrt(4s²/(ad(1+as²)² - (1-as²)²))
    x = (s + s) * Dx
    x_neg = is_negative(x)
    x = conditional_negate(x, x_neg)

    # y == (1-as²)/(1+as²)
    y = u1 * Dy

    # t == ((1+as²) sqrt(4s²/(ad(1+as²)² - (1-as²)²)))/(1-as²)
    t = x * y

    if !ok || is_negative(t) || iszero(y)
        nothing
    else
        RistrettoPoint(EdwardsPoint(x, y, one_, t))
    end
end


struct RistrettoPoint{T}
    ep :: EdwardsPoint{T}
end


Base.zero(::Type{RistrettoPoint{T}}) where T = RistrettoPoint{T}(zero(EdwardsPoint{T}))


# Compress this point using the Ristretto encoding.
function compress(p::EdwardsPoint{T}) where T
    X = p.ep.X
    Y = p.ep.Y
    Z = p.ep.Z
    T_ = p.ep.T_

    u1 = (Z + Y) * (Z - Y)
    u2 = X * Y
    # Ignore return value since this is always square
    (_, invsqrt_) = invsqrt(u1 * square(u2))
    i1 = invsqrt_ * u1
    i2 = invsqrt_ * u2
    z_inv = i1 * (i2 * T_)
    den_inv = i2

    iX = X * SQRT_M1
    iY = Y * SQRT_M1
    ristretto_magic = INVSQRT_A_MINUS_D
    enchanted_denominator = i1 * ristretto_magic

    rotate = is_negative(T_ * z_inv)

    X = conditional_assign(X, iY, rotate)
    Y = conditional_assign(Y, iX, rotate)
    den_inv = conditional_assign(den_inv, enchanted_denominator, rotate)

    Y = conditional_negate(Y, is_negative(X * z_inv))

    s = den_inv * (Z - Y)
    s_is_negative = is_negative(s)
    s = conditional_negate(s, s_is_negative)

    CompressedRistretto{T}(s)
end


Base.:+(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T = RistrettoPoint{T}(p.ep + q.ep)


Base.:*(p::RistrettoPoint{T}, s::Z) where {T, Z} = RistrettoPoint{T}(p.ep * s)


function Base.:(==)(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T
    X1Y2 = p.ep.X * p.ep.Y
    Y1X2 = p.ep.Y * p.ep.X
    X1X2 = p.ep.X * p.ep.X
    Y1Y2 = p.ep.Y * p.ep.Y

    ct_eq(X1Y2, Y1X2) | ct_eq(X1X2, Y1Y2)
end
