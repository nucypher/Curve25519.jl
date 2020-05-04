#=
Ported from backend/serial/curve_models/mod.rs
=#


@inline square(x) = x * x


# Returns 2 times the square of this field element.
@inline square2(x) = let s = square(x)
    s + s
end


struct ProjectivePoint{T}
    X :: T
    Y :: T
    Z :: T
end


@inline Base.zero(::Type{ProjectivePoint{T}}) where T = ProjectivePoint{T}(zero(T), one(T), one(T))


# Double this point: return self + self
@inline function double(p::ProjectivePoint{T}) where T
    XX          = square(p.X)
    YY          = square(p.Y)
    ZZ2         = square2(p.Z)
    X_plus_Y    = p.X + p.Y
    X_plus_Y_sq = square(X_plus_Y)
    YY_plus_XX  = YY + XX
    YY_minus_XX = YY - XX

    CompletedPoint{T}(X_plus_Y_sq - YY_plus_XX, YY_plus_XX, YY_minus_XX, ZZ2 - YY_minus_XX)
end


#=
Convert this point from the \\( \mathbb P\^2 \\) model to the
\\( \mathbb P\^3 \\) model.

This costs \\(3 \mathrm M + 1 \mathrm S\\).
=#
@inline function to_extended(p::ProjectivePoint{T}) where T
    EdwardsPoint{T}(p.X * p.Z, p.Y * p.Z, square(p.Z), p.X * p.Y)
end


#=
A `CompletedPoint` is a point \\(((X:Z), (Y:T))\\) on the
\\(\mathbb P\^1 \times \mathbb P\^1 \\) model of the curve.
A point (x,y) in the affine model corresponds to \\( ((x:1),(y:1)) \\).
=#
struct CompletedPoint{T}
    X :: T
    Y :: T
    Z :: T
    T_ :: T
end


#=
Convert this point from the
\\( \mathbb P\^1 \times \mathbb P\^1 \\) model to the \\( \mathbb P\^2 \\) model.

This costs \\(3 \mathrm M \\).
=#
@inline function to_projective(p::CompletedPoint{T}) where T
    ProjectivePoint(p.X * p.T_, p.Y * p.Z, p.Z * p.T_)
end


#=
Convert this point from the
\\( \mathbb P\^1 \times \mathbb P\^1 \\) model to the \\( \mathbb P\^3 \\) model.

This costs \\(4 \mathrm M \\).
=#
@inline function to_extended(p::CompletedPoint{T}) where T
    EdwardsPoint{T}(p.X * p.T_, p.Y * p.Z, p.Z * p.T_, p.X * p.Y)
end


#=
A pre-computed point on the \\( \mathbb P\^3 \\) model for the
curve, represented as \\((Y+X, Y-X, Z, 2dXY)\\) in "Niels coordinates".
=#
struct ProjectiveNielsPoint{T} <: CT.Selectable
    Y_plus_X :: T
    Y_minus_X :: T
    Z :: T
    T2d :: T
end


@inline function CT.select(
        choice::CT.Choice, p::ProjectiveNielsPoint{T}, q::ProjectiveNielsPoint{T}) where T <: CT.Value
    ProjectiveNielsPoint{T}(
        CT.select(choice, p.Y_plus_X, q.Y_plus_X),
        CT.select(choice, p.Y_minus_X, q.Y_minus_X),
        CT.select(choice, p.Z, q.Z),
        CT.select(choice, p.T2d, q.T2d)
        )
end


@inline Base.zero(::Type{ProjectiveNielsPoint{T}}) where T = ProjectiveNielsPoint{T}(one(T), one(T), one(T), zero(T))


@inline function Base.:-(p::ProjectiveNielsPoint{T}) where T
    ProjectiveNielsPoint{T}(p.Y_minus_X, p.Y_plus_X, p.Z, -p.T2d)
end


@inline function Base.:+(self::EdwardsPoint{T}, other::ProjectiveNielsPoint{T}) where T
    Y_plus_X  = self.Y + self.X
    Y_minus_X = self.Y - self.X
    PP = Y_plus_X  * other.Y_plus_X
    MM = Y_minus_X * other.Y_minus_X
    TT2d = self.T_ * other.T2d
    ZZ   = self.Z * other.Z
    ZZ2  = ZZ + ZZ

    CompletedPoint{T}(PP - MM, PP + MM, ZZ2 + TT2d, ZZ2 - TT2d)
end


struct AffineNielsPoint{T}
    y_plus_x :: T
    y_minus_x :: T
    xy2d :: T
end


@inline Base.zero(::Type{AffineNielsPoint{T}}) where T = AffineNielsPoint{T}(one(T), one(T), zero(T))


function Base.:+(p::EdwardsPoint{T}, q::AffineNielsPoint{T}) where T
    Y_plus_X  = p.Y + p.X
    Y_minus_X = p.Y - p.X
    PP        = Y_plus_X  * q.y_plus_x
    MM        = Y_minus_X * q.y_minus_x
    Txy2d     = p.T_ * q.xy2d
    Z2        = p.Z + p.Z

    CompletedPoint{T}(PP - MM, PP + MM, Z2 + Txy2d, Z2 - Txy2d)
end


@inline function Base.:-(p::AffineNielsPoint{T}) where T
    AffineNielsPoint{T}(p.y_minus_x, p.y_plus_x, -p.xy2d)
end
