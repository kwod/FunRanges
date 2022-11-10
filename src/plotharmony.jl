using Plots, FunRanges
fr = FunRange(1/5,5,800,logarithmic)
harmony(r::Rational) = inv(√(abs(r.num) * r.den))
nds = bestrats(fr)
rs = Rational.(nds)
pars = @. Rational(parentnode(nds))
p = sticks(rs, harmony.(rs), xscale = :log10, legend = false)
for (p1, p2) in zip(rs, pars)
    if .2 <= p2 <= 5
        color = p1 < p2 ? :red : :green
        plot!(p, [(p1, harmony(p1)),(p2, harmony(p2))], color = color)
    end
end
p