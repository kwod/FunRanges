using FunRanges
using Test

@testset "FunRanges.jl" begin
    fr = FunRange(1, 2, 13, logarithmic)
    @test length(fr) == 13
end
