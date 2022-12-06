# SPDX-License-Identifier: MIT

struct SSBTNode
    lo::Vector{BigInt}
    hi::Vector{BigInt}
    lev::Int
    num::BigInt
end

function mid(nd::SSBTNode)
    iszero(nd.lev) ? [big(0), big(1)] : nd.lo .+ nd.hi
end

function Base.Rational(nd::SSBTNode)
    Rational(mid(nd)...)
end

SSBT_root = SSBTNode(big.([-1, 0]), big.([1, 0]), 0, big(0))

function parentnode(nd::SSBTNode)
    iszero(nd.lev) ? nd :
        isodd(nd.num) ?
            SSBTNode(nd.lo .- nd.hi, nd.hi, nd.lev - 1, nd.num >> 1) :
            SSBTNode(nd.lo, nd.hi .- nd.lo, nd.lev - 1, nd.num >> 1)
end

function childnode(nd::SSBTNode, right::Bool)
    lo, hi = right ? (mid(nd), nd.hi) : (nd.lo, mid(nd))
    SSBTNode(lo, hi, nd.lev + 1, 2 * nd.num + Int(right))
end

function SSBTNode(x::Union{Rational, Integer, AbstractFloat})
    if isfinite(x)
        reduce(childnode, bits(x), init = SSBT_root)
    end
end

function id(nd::SSBTNode)
    big(2)^nd.lev + nd.num
end

function id(x::Union{Rational, Integer, AbstractFloat})
    id(SSBTNode(x))
end

function qm(nd::SSBTNode)
    (nd.num + 0.5) / big(2)^nd.lev
end

function qm(x::Union{Rational, Integer, AbstractFloat})
    qm(SSBTNode(x))
end