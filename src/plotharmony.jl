using Plots, FunRanges
fr1 = FunRange(1,2,1300,logarithmic)
harmony(r::Rational) = inv(√(abs(r.num) * r.den))
br1, pr1 = bestrats(fr1)
p = sticks(br1, harmony.(br1), xscale = :log10, legend = false)
for (p1, p2) in zip(br1, pr1)
    if 1 <= p2 <= 2
        color = p1 < p2 ? :red : :green
        plot!(p, [(p1, harmony(p1)),(p2, harmony(p2))], color = color)
    end
end
p