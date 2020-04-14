# in Dalek it is just `ProjectivePoint`, clashing from the one from montgomery.jl
# Hopefully I can distinguish which is which by the number of arguments
struct IntProjectivePoint{T}
    X :: T
    Y :: T
    Z :: T
end


Base.zero(::Type{IntProjectivePoint{T}}) where T = IntProjectivePoint{T}(zero(T), one(T), one(T))


# Double this point: return self + self
function double(p::IntProjectivePoint{T}) where T
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
/// Convert this point from the \\( \mathbb P\^2 \\) model to the
/// \\( \mathbb P\^3 \\) model.
///
/// This costs \\(3 \mathrm M + 1 \mathrm S\\).
=#
function to_extended(p::IntProjectivePoint{T}) where T
    EdwardsPoint{T}(p.X * p.Z, p.Y * p.Z, square(p.Z), p.X * p.Y)
end


#=
/// A `CompletedPoint` is a point \\(((X:Z), (Y:T))\\) on the \\(\mathbb
/// P\^1 \times \mathbb P\^1 \\) model of the curve.
/// A point (x,y) in the affine model corresponds to \\( ((x:1),(y:1))
/// \\).
///
/// More details on the relationships between the different curve models
/// can be found in the module-level documentation.
#[derive(Copy, Clone)]
#[allow(missing_docs)]
=#
struct CompletedPoint{T}
    X :: T
    Y :: T
    Z :: T
    T_ :: T
end


#=
/// Convert this point from the \\( \mathbb P\^1 \times \mathbb P\^1
/// \\) model to the \\( \mathbb P\^2 \\) model.
///
/// This costs \\(3 \mathrm M \\).
=#
function to_projective(p::CompletedPoint{T}) where T
    IntProjectivePoint(p.X * p.T_, p.Y * p.Z, p.Z * p.T_)
end


#=
/// Convert this point from the \\( \mathbb P\^1 \times \mathbb P\^1
/// \\) model to the \\( \mathbb P\^3 \\) model.
///
/// This costs \\(4 \mathrm M \\).
=#
function to_extended(p::CompletedPoint{T}) where T
    EdwardsPoint{T}(p.X * p.T_, p.Y * p.Z, p.Z * p.T_, p.X * p.Y)
end


#=
/// A pre-computed point on the \\( \mathbb P\^3 \\) model for the
/// curve, represented as \\((Y+X, Y-X, Z, 2dXY)\\) in "Niels coordinates".
///
/// More details on the relationships between the different curve models
/// can be found in the module-level documentation.
#[derive(Copy, Clone)]
=#
struct ProjectiveNielsPoint{T}
    Y_plus_X :: T
    Y_minus_X :: T
    Z :: T
    T2d :: T
end


Base.zero(::Type{ProjectiveNielsPoint{T}}) where T = ProjectiveNielsPoint{T}(one(T), one(T), one(T), zero(T))


function Base.:-(p::ProjectiveNielsPoint{T}) where T
    ProjectiveNielsPoint{T}(p.Y_minus_X, p.Y_plus_X, p.Z, -p.T2d)
end


function Base.:+(self::EdwardsPoint{T}, other::ProjectiveNielsPoint{T}) where T
    Y_plus_X  = self.Y + self.X
    Y_minus_X = self.Y - self.X
    PP = Y_plus_X  * other.Y_plus_X
    MM = Y_minus_X * other.Y_minus_X
    TT2d = self.T_ * other.T2d
    ZZ   = self.Z * other.Z
    ZZ2  = ZZ + ZZ

    CompletedPoint{T}(PP - MM, PP + MM, ZZ2 + TT2d, ZZ2 - TT2d)
end

