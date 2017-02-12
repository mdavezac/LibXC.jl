module LibXCTests
using LibXC
using Base.Test

@testset "> Internal API" begin
    @testset ">> xc key from name" begin
        @test LibXC.version == v"3.0.0"
    end
end
end
