module LibXCTests
using LibXC
using Base.Test

@testset "> Internal API" begin
    @test LibXC.LIB_VERSION == v"3.0.0"
    @test LibXC.FUNCTIONALS[:lda_x] == 1
    @test LibXC.FUNCTIONALS[:mgga_x_ms2] == 223
end


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
@testset "> Internal Functional API: $symb" for (symb, stuff) in functionals
    func = LibXC.XCFunctional(symb, stuff[:spin])
    info_ptr = unsafe_load(unsafe_load(func.c_ptr).info)
    @test unsafe_string(info_ptr.name) == stuff[:description]
    @test description(func) == stuff[:description]
    @test family(func) == stuff[:family]
    @test kind(func) == stuff[:kind]
    @test LibXC.libkey(func) == stuff[:libkey]
    @test LibXC.spin(func) == stuff[:spin]
end


end
