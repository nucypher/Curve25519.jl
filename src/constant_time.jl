#=
Stubs for constant-time conditional operations.
TODO: needs to be constant-time.
=#

ct_eq(x, y) = x == y


conditional_assign(x, y, choice) = choice ? y : x


conditional_negate(x, choice) = choice ? -x : x
