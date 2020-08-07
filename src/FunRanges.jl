module FunRanges
import Base: iterate, length, getindex, show, convert, first, last
export FunType, FunRange, bestrats, logarithmic, circular

@enum FunType linear logarithmic circular

fun = Dict(linear => identity, logarithmic => log, circular => atan)
rev = Dict(linear => identity, logarithmic => exp, circular => tan)

struct FunRange{T} <: AbstractRange{T}
    lr::LinRange{T}
    ft::FunType
end

function FunRange(start, stop, len, ft)
    FunRange(
        LinRange(fun[ft](start), fun[ft](stop), len),
        ft)
end

function iterate(r::FunRange, state = iterate(r.lr))
    state == nothing && return nothing
    val, next = state
    rev[r.ft](val), iterate(r.lr, next)
end

length(r::FunRange) = length(r.lr)

getindex(r::FunRange, s::OrdinalRange) = rev[r.ft].(getindex(r.lr, s))
getindex(r::FunRange, i::Integer) = rev[r.ft].(getindex(r.lr, i))

first(r::FunRange) = rev[r.ft](first(r.lr))

last(r::FunRange) = rev[r.ft](last(r.lr))

show(io::IO, r::FunRange) = print(io, "$(first(r)):$(length(r)-1)_$(r.ft)_steps:$(last(r))")

function findindex(r::FunRange, val::Real)
    q = (fun[r.ft](val) - first(r.lr)) / step(r.lr) + 1
    q = max(q, 0)
    q = min(q, length(r)+1)
    round(Int, q)
end

convert(newft::FunType, r::FunRange) = FunRange(first(r), last(r), length(r), newft)

function zoom(r::FunRange, zoomfactor)
    delta = (1-inv(zoomfactor)) * (last(r.lr) - first(r.lr)) / 2
    FunRange(LinRange(r.lr.start+delta, r.lr.stop-delta, length(r)), r.ft)
end

function shift(r::FunRange, shiftpart)
    delta = shiftpart * (last(r.lr) - first(r.lr))
    FunRange(LinRange(r.lr.start+delta, r.lr.stop+delta, length(r)), r.ft)
end

function bestrats(fr::FunRange)
    rats = similar(Array{Rational{Int}}, (length(fr),))
    function approxrats(lo, hi, i, j)
        mid = lo .+ hi
        mid == [0, 0] && (mid = [0, 1])
        rat = Rational(mid...)
        k = findindex(fr, rat)
        if k > j
            approxrats(lo, mid, i, j)
        elseif k < i
            approxrats(mid, hi, i, j)
        else
            @inbounds rats[k] = rat
            k > i && approxrats(lo, mid, i, k - 1)
            k < j && approxrats(mid, hi, k + 1, j)
        end
    end
    approxrats([-1, 0], [1, 0], 1, length(fr))
    rats
end

end # module