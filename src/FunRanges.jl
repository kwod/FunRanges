module FunRanges

include("Bits.jl")
include("SSBTNodes.jl")
include("harmony.jl")

import Base: iterate, length, getindex, show, convert, first, last, Rational
using Base.Threads, Plots
export FunType, FunRange, bestrats, linear, logarithmic, circular,
    SSBTNode, childnode, parentnode, qm, id, bits,
    harmony, plotharmony

@enum FunType linear logarithmic circular

fun = Dict(linear => identity, logarithmic => log, circular => atan)
rev = Dict(linear => identity, logarithmic => exp, circular => tan)

struct FunRange{T} <: AbstractRange{T}
    lr::LinRange{T}
    ft::FunType
end

function FunRange(start::Real, stop::Real, len::Real, ft::FunType)
    FunRange(
        LinRange(fun[ft](start), fun[ft](stop), len),
        ft)
end

function iterate(r::FunRange, state = iterate(r.lr))
    if !isnothing(state)
        val, next = state
        rev[r.ft](val), iterate(r.lr, next)
    end
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

function zoom(r::FunRange, zoomfactor::Real)
    delta = (1-inv(zoomfactor)) * (last(r.lr) - first(r.lr)) / 2
    FunRange(LinRange(r.lr.start+delta, r.lr.stop-delta, length(r)), r.ft)
end

function shift(r::FunRange, shiftpart::Real)
    delta = shiftpart * (last(r.lr) - first(r.lr))
    FunRange(LinRange(r.lr.start+delta, r.lr.stop+delta, length(r)), r.ft)
end

function bestrats(fr::FunRange)
    nds = similar(SSBTNode[], length(fr))
    ancestorindex = similar(Int[], length(fr))
    function approxrats(nd, i, j, anci)
        rat = Rational(nd)
        k = findindex(fr, rat)
        if k > j
            approxrats(childnode(nd, false), i, j, anci)
        elseif k < i
            approxrats(childnode(nd, true), i, j, anci)
        else
            @inbounds nds[k] = nd
            @inbounds ancestorindex[k] = anci
            if k > i 
                lt = @spawn approxrats(childnode(nd, false), i, k - 1, k)
            end
            k < j && approxrats(childnode(nd, true), k + 1, j, k)
            k > i && wait(lt)      
        end
    end
    approxrats(SSBT_root, 1, length(fr), 0)
    Rational.(nds), ancestorindex
end

end # module