functionals =  [
    :lda_x => Dict{Symbol, Any}(
        :description=>"Slater exchange",
        :family => LibXC.Constants.lda,
        :kind => LibXC.Constants.exchange,
        :libkey => 1,
        :flags => Set(LibXC.Constants.FLAGS[
            LibXC.Constants.kxc, LibXC.Constants.exc, LibXC.Constants.D3,
            LibXC.Constants.vxc, LibXC.Constants.fxc]),
        :spin => convert(LibXC.Constants.SPIN, 1)
    ), :hyb_mgga_xc_m11 => Dict{Symbol, Any}(
        :description => "Minnesota M11 hybrid functional",
        :family => LibXC.Constants.hybrid_mgga,
        :kind => LibXC.Constants.excorr,
        :libkey => 462,
        :flags => Set(LibXC.Constants.FLAGS[
            LibXC.Constants.kxc, LibXC.Constants.exc, LibXC.Constants.D3,
            LibXC.Constants.vxc, LibXC.Constants.fxc]),
        :spin => convert(LibXC.Constants.SPIN, 2)
    )
]
@testset ">> $symb" for (symb, stuff) in functionals
    func = LibXC.XCFunctional(symb, stuff[:spin])
    info_ptr = unsafe_load(unsafe_load(func.c_ptr).info)
    @test unsafe_string(info_ptr.name) == stuff[:description]
    @test description(func) == stuff[:description]
    @test family(func) == stuff[:family]
    @test kind(func) == stuff[:kind]
    @test LibXC.Internals.libkey(func) == stuff[:libkey]
    @test LibXC.spin(func) == stuff[:spin]

#     @test LibXC.output_size(LibXC.Constants.polarized, (2, 5), 1) == (5,)
#     @test LibXC.output_size(LibXC.Constants.unpolarized, (2, 5), 1) == (2, 5)
#     @test LibXC.output_size(LibXC.Constants.polarized, (6,), 1) == (3,)
#     @test LibXC.output_size(LibXC.Constants.unpolarized, (6,), 1) == (6,)
#     @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (5,), 1)
#     @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (), 1)
#     @test_throws ArgumentError LibXC.output_size(LibXC.Constants.unpolarized, (), 1)
#     @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (3, 2), 1)
end
