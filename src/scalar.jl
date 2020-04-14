bot_half(x::UInt8) = (x >> 0) & 0xf
top_half(x::UInt8) = (x >> 4) & 0xf


#=
/// Write this scalar in radix 16, with coefficients in \\([-8,8)\\),
/// i.e., compute \\(a\_i\\) such that
/// $$
///    a = a\_0 + a\_1 16\^1 + \cdots + a_{63} 16\^{63},
/// $$
/// with \\(-8 \leq a_i < 8\\) for \\(0 \leq i < 63\\) and \\(-8 \leq a_{63} \leq 8\\).
=#
function to_radix_16(x::T) where T
    # debug_assert!(self[31] <= 127); - we don't need this if T is already modulo 2^255-19
    output = zeros(Int8, 64)

    # TODO: can be sped up
    ubytes = [convert(UInt8, (x >> ((i-1) * 8)) & 0xff) for i in 1:32]

    # Step 1: change radix.
    # Convert from radix 256 (bytes) to radix 16 (nibbles)
    for i in 0:31
        output[2*i+1  ] = signed(bot_half(ubytes[i+1]))
        output[2*i+1+1] = signed(top_half(ubytes[i+1]))
    end
    # Precondition note: since self[31] <= 127, output[63] <= 7

    # Step 2: recenter coefficients from [0,16) to [-8,8)
    for i in 0:62
        carry    = (output[i+1] + 8) >> 4
        output[i+1  ] -= carry << 4
        output[i+1+1] += carry
    end
    # Precondition note: output[63] is not recentered.  It
    # increases by carry <= 1.  Thus output[63] <= 8.

    output
end


function to_radix_16(x::AbstractModUInt)
    to_radix_16(value(x))
end
