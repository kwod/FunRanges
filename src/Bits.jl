struct Bits{T}
    a::T
    b::T
end

function Base.iterate(r::Bits, (a, b) = (r.a, r.b))
    if a < b
        false, (a, b - a)
    elseif a > b
        true, (a - b, b)
    end
end

function Base.eltype(::Type{Bits})
    Bool
end

function Base.IteratorSize(::Type{Bits})
    Base.SizeUnknown()
end

"""
    bits(r::Rational)

bijectiv mapping of real numbers to iterated bits
"""
function bits(r::Rational{T}) where T
    if isfinite(r)
        Iterators.flatten(
            if r < 0
                (false, Bits(r.den, -r.num))
            elseif r > 0
                (true, Bits(r.num, r.den))
            else
                ()
            end
            )
    end
end

function bits(i::Integer)
    bits(i//one(i))
end

function bits(flt::AbstractFloat)
    bits(rationalize(flt))
end