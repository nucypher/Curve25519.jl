struct CompressedEdwardsY{T}
    y :: T

    CompressedEdwardsY{T}(x::T) where T = new{T}(x)
    CompressedEdwardsY(x::T) where T = new{T}(x)
end


function decompress(p::CompressedEdwardsY{T}) where T
    Y = p.Y
    Z = one(T)
    YY = square(Y)
    u = YY - Z # u = y²-1
    v = (YY * EDWARDS_D) + Z # v = dy²+1
    (is_valid_y_coord, X) = sqrt_ratio_i(u, v)

    if !is_valid_y_coord
        return nothing
    end

    # FieldElement::sqrt_ratio_i always returns the nonnegative square root,
    # so we negate according to the supplied sign bit.
    compressed_sign_bit = !iszero(p.y >> 255) # Choice::from(self.as_bytes()[31] >> 7);
    X = conditional_negate(X, compressed_sign_bit)

    EdwardsPoint{T}(X, Y, Z, X * Y)
end


struct EdwardsPoint{T}
    X :: T
    Y :: T
    Z :: T
    T_ :: T

    EdwardsPoint{T}(x::T, y::T, z::T, t::T) where T = new{T}(x, y, z, t)
    EdwardsPoint(x::T, y::T, z::T, t::T) where T = new{T}(x, y, z, t)
end


Base.zero(::Type{EdwardsPoint{T}}) where T = EdwardsPoint{T}(zero(T), one(T), one(T), zero(T))


function is_valid(p::EdwardsPoint{T}) where T
    point_on_curve = is_valid(to_projective(p))
    on_segre_image = (p.X * p.Y) == (p.Z * p.T_)

    point_on_curve && on_segre_image
end


function to_projective_niels(p::EdwardsPoint{T}) where T
    ProjectiveNielsPoint{T}(p.Y + p.X, p.Y - p.X, p.Z, p.T_ * EDWARDS_D2)
end


function to_projective(p::EdwardsPoint{T}) where T
    IntProjectivePoint{T}(p.X, p.Y, p.Z)
end


double(p::EdwardsPoint{T}) where T = to_extended(double(to_projective(p)))


Base.:+(p::EdwardsPoint{T}, q::EdwardsPoint{T}) where T = to_extended(p + to_projective_niels(q))


Base.:-(p::EdwardsPoint{T}) where T = EdwardsPoint{T}(-p.X, p.Y, p.Z, -p.T)


function Base.:*(point::EdwardsPoint{T}, scalar::Z) where {T, Z}
    scalar_mul__variable_base__mul(point, scalar)
end

