# TODO: needs to be constant-time. See ct_eq() implementation in Dalek

ct_eq(x, y) = x == y

conditional_assign(x, y, choice) = choice ? y : x

conditional_negate(x, choice) = choice ? -x : x

function conditional_swap(p, q, choice)
    choice ? (q, p) : (p, q)
end


struct MontgomeryPoint{T}
    x :: T
end


# in Dakek: default()
Base.zero(::Type{MontgomeryPoint{T}}) where T = MontgomeryPoint{T}(zero(T))


Base.:(==)(p::MontgomeryPoint{T}, q::MontgomeryPoint{T}) where T = ct_eq(p.x, q.x)


function to_edwards(p::MontgomeryPoint{T}, sign_::Bool) where T
    #=
    // To decompress the Montgomery u coordinate to an
    // `EdwardsPoint`, we apply the birational map to obtain the
    // Edwards y coordinate, then do Edwards decompression.
    //
    // The birational map is y = (u-1)/(u+1).
    //
    // The exceptional points are the zeros of the denominator,
    // i.e., u = -1.
    //
    // But when u = -1, v^2 = u*(u^2+486662*u+1) = 486660.
    //
    // Since this is nonsquare mod p, u = -1 corresponds to a point
    // on the twist, not the curve, so we can reject it early.
    =#
    u = p.x

    if u == minus_one(T)
        return nothing
    end

    one_ = one(T)

    y = (u - one_) * inv(u + one_)
    if sign_
        y ^= one(T) << 255 # TODO: check that it's equivalent to `y_bytes[31] ^= sign << 7;`
    end

    decompress(CompressedEdwardsY(y))
end


struct ProjectivePoint{T}
    U :: T
    W :: T
end


# In Dalek: identity() and default()
Base.zero(::Type{ProjectivePoint{T}}) where T = ProjectivePoint{T}(one(T), zero(T))


function to_affine(p::ProjectivePoint{T}) where T
    MontgomeryPoint{T}(p.U * inv(p.W))
end


square(x) = x * x


function differential_add_and_double(p::P, q::P, affine_PmQ::T) where {P <: ProjectivePoint{T}} where T
    t0 = p.U + p.W
    t1 = p.U - p.W
    t2 = q.U + q.W
    t3 = q.U - q.W

    t4 = square(t0) # (U_P + W_P)^2 = U_P^2 + 2 U_P W_P + W_P^2
    t5 = square(t1) # (U_P - W_P)^2 = U_P^2 - 2 U_P W_P + W_P^2

    t6 = t4 - t5 # 4 U_P W_P

    t7 = t0 * t3 # (U_P + W_P) (U_Q - W_Q) = U_P U_Q + W_P U_Q - U_P W_Q - W_P W_Q
    t8 = t1 * t2 # (U_P - W_P) (U_Q + W_Q) = U_P U_Q - W_P U_Q + U_P W_Q - W_P W_Q

    t9 = t7 + t8 # 2 (U_P U_Q - W_P W_Q)
    t10 = t7 - t8 # 2 (W_P U_Q - U_P W_Q)

    t11 = square(t9) # 4 (U_P U_Q - W_P W_Q)^2
    t12 = square(t10) # 4 (W_P U_Q - U_P W_Q)^2

    t13 = APLUS2_OVER_FOUR * t6 # (A + 2) U_P U_Q

    t14 = t4 * t5 # ((U_P + W_P)(U_P - W_P))^2 = (U_P^2 - W_P^2)^2
    t15 = t13 + t5 # (U_P - W_P)^2 + (A + 2) U_P W_P

    t16 = t6 * t15 # 4 (U_P W_P) ((U_P - W_P)^2 + (A + 2) U_P W_P)

    t17 = affine_PmQ * t12 # U_D * 4 (W_P U_Q - U_P W_Q)^2
    t18 = t11 # W_D * 4 (U_P U_Q - W_P W_Q)^2

    #=
    U_{P'} = (U_P + W_P)^2 (U_P - W_P)^2
    W_{P'} = (4 U_P W_P) ((U_P - W_P)^2 + ((A + 2)/4) 4 U_P W_P)
    U_{Q'} = W_D * 4 (U_P U_Q - W_P W_Q)^2
    W_{Q'} = U_D * 4 (W_P U_Q - U_P W_Q)^2
    =#
    ProjectivePoint{T}(t14, t16), ProjectivePoint{T}(t18, t17)
end


function Base.:*(p::MontgomeryPoint{T}, s::Integer) where T
    # Algorithm 8 of Costello-Smith 2017

    x0 = zero(ProjectivePoint{T})
    x1 = ProjectivePoint{T}(p.X, one(T))

    # A slight devication from Costello-Smith to avoid defining and using
    # a separate double() function to initialize x0.
    nb = num_bits(s) + 1

    bits = [isodd(s >> i) for i in 0:nb-1]

    for i in nb-2:-1:0
        choice = xor(bits[i+2], bits[i+1])
        x0, x1 = conditional_swap(x0, x1, choice)
        x0, x1 = add_and_double(x0, x1, p.X)
    end
    x0, x1 = conditional_swap(x0, x1, bits[1])
    to_affine(x0)
end


function is_valid(p::ProjectivePoint{T}) where T
    # Curve equation is    -x^2 + y^2 = 1 + d*x^2*y^2,
    # homogenized as (-X^2 + Y^2)*Z^2 = Z^4 + d*X^2*Y^2
    XX = square(p.X)
    YY = square(p.Y)
    ZZ = square(p.Z)
    ZZZZ = square(ZZ)
    lhs = (YY - XX) * ZZ
    rhs = ZZZZ + EDWARDS_D * (XX * YY)

    lhs == rhs
end
