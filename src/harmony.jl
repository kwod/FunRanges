harmony(r::Rational{T}) where T = inv(√(abs(r.num) * r.den))

function plotharmony()
    lns = 7000
    fr = FunRange(.5,4,lns,logarithmic)
    rs, ps = bestrats(fr)
    xs, ys = Float64.(rs), harmony.(rs)
    p = plot(xscale = :log10, legend = false)
    for i in 1:lns
        j = ps[i]
        if !iszero(j)
            color = i < j ? :red : :green
            plot!(p, [xs[i], xs[j]], [ys[i], ys[j]], color = color)
        end
    end
    sticks!(p, xs, ys; color = :blue)
    p
end