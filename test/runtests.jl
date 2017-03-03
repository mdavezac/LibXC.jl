module LibXCTests
using LibXC
using Base.Test
using CSV: read
using DataFrames: DataFrame

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

    @test LibXC.εxc_size(true, (2, 5)) == (5,)
    @test LibXC.εxc_size(false, (2, 5)) == (2, 5)
    @test LibXC.εxc_size(true, (6,)) == (3,)
    @test LibXC.εxc_size(false, (6,)) == (6,)
    @test_throws ArgumentError LibXC.εxc_size(true, (5,))
    @test_throws ArgumentError LibXC.εxc_size(true, ())
    @test_throws ArgumentError LibXC.εxc_size(false, ())
    @test_throws ArgumentError LibXC.εxc_size(true, (3, 2))
end

files = ["lda_x.BrOH.unpol.dat", "lda_x.BrOH.pol.dat"]
@testset "> LDA $name" for name in files
    input_data_file = joinpath(dirname(Base.source_path()), "data", split(name, '.')[2])
    expected_data_file = joinpath(dirname(Base.source_path()), "data", name)

    input = (map(x -> parse(Float64, x), split(v))
             for v in readlines(input_data_file)[2:end])
    input = transpose(hcat(input...))
    input = DataFrame(Any[input[:, i] for i in 1:size(input, 2)],
                      [:ρ_a, :ρ_b, :σ_aa, :σ_ab, :σ_bb, :δ_a, :δ_b, :τ_a, :τ_b])

    expected = (map(x -> parse(Float64, x), split(v))
               for v in readlines(expected_data_file)[3:end])
    expected = transpose(hcat(expected...))
    if size(expected, 2) == 3
        expected = DataFrame(Any[expected[:, i] for i in 1:size(expected, 2)],
                             [:ε, :v, :δv])

        ρ = input[:ρ_a] + input[:ρ_b]
        @test energy(:lda_x, ρ) ≈ expected[:ε]
    else
        expected = DataFrame(Any[expected[:, i] for i in 1:size(expected, 2)],
                             [:ε, :v_a, :v_b, :δv_aa, :δv_ab, :δv_bb])

        ρs = vcat(input[:ρ_a]', input[:ρ_b]')
        @test energy(:lda_x, ρs) ≈ expected[:ε]
    end
end



end
