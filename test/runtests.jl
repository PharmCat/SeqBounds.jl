using SeqBounds
using Test

@testset "SeqBounds.jl" begin
    io = IOBuffer();
    b = SeqBounds.bounds([0.2, 0.4, 0.6, 0.8, 1.0], 0.05; h = 0.05, side = :one, asf = :obf)
    @test b.zb[end] ≈ 1.7396609804214789 atol=1E-6
    @test_nowarn Base.show(io, b)
    b = SeqBounds.bounds([0.2, 0.4, 0.6, 0.8, 1.0], 0.05; h = 0.05, side = :two, asf = :pocock)
    @test b.zb[end] ≈ 2.385909219073807 atol=1E-6
    @test_nowarn Base.show(io, b)

    @test_nowarn SeqBounds.bounds([0.25, 0.5, 0.75, 1.0], 0.05; h = 0.04, side = :two, asf = :power)
    @test_nowarn SeqBounds.bounds([0.25, 0.5, 0.75, 1.0], 0.025; h = 0.04, side = :two, asf = :ep)
end
