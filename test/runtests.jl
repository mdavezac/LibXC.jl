module LibXCTests
using LibXC
using Base.Test

@testset "> Internal API" begin
    @testset ">> xc key from name" begin
        @test LibXC.lib_version == v"3.0.0"
        @test LibXC._xc_key("LDA_X") == 1
        @test LibXC._xc_key("lda_x") == 1
        @test LibXC._xc_name(1) == "lda_x"
    end
end
end
