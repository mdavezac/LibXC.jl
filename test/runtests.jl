module LibXCTests
using LibXC
using Base.Test

@testset "> Internal API" begin
    @test LibXC.LIB_VERSION == v"3.0.0"
    @test LibXC.FUNCTIONALS[:lda_x] == 1
    @test LibXC.FUNCTIONALS[:mgga_x_ms2] == 223
end


functionals =  [:lda_x => "Slater exchange",
                :hyb_mgga_xc_m11 => "Minnesota M11 hybrid functional"]
@testset "> Internal Functional API: $symb" for (symb, name) in functionals
    func = LibXC.XCFunctional(symb, true)
    info_ptr = unsafe_load(unsafe_load(func.c_ptr).info)
    @test unsafe_string(info_ptr.name) == name
end

end
