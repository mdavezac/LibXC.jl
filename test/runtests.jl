module LibXCTests
using LibXC
using Base.Test
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

    @test LibXC.output_size(LibXC.Constants.polarized, (2, 5), 1) == (5,)
    @test LibXC.output_size(LibXC.Constants.unpolarized, (2, 5), 1) == (2, 5)
    @test LibXC.output_size(LibXC.Constants.polarized, (6,), 1) == (3,)
    @test LibXC.output_size(LibXC.Constants.unpolarized, (6,), 1) == (6,)
    @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (5,), 1)
    @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (), 1)
    @test_throws ArgumentError LibXC.output_size(LibXC.Constants.unpolarized, (), 1)
    @test_throws ArgumentError LibXC.output_size(LibXC.Constants.polarized, (3, 2), 1)
end

function input_data(name::String)
    data_file = joinpath(dirname(Base.source_path()), "data", name)

    input = (map(x -> parse(Float64, x), split(v))
             for v in readlines(data_file)[2:end])
    input = transpose(hcat(input...))
    input = DataFrame(Any[input[:, i] for i in 1:size(input, 2)],
                      [:ρ_a, :ρ_b, :σ_aa, :σ_ab, :σ_bb, :δ_a, :δ_b, :τ_a, :τ_b])
end

function expected_data(name::String)
    data_file = joinpath(dirname(Base.source_path()), "data", name)
    expected = (map(x -> parse(Float64, x), split(v))
                for v in readlines(data_file)[3:end])
    expected = transpose(hcat(expected...))
    if size(expected, 2) == 3
        DataFrame(Any[reinterpret(LibXC.EnergyDensity{Cdouble}, expected[:, i])
                      for i in 1:size(expected, 2)], [:ε, :v, :δv])
    elseif size(expected, 2) == 6
        DataFrame(Any[reinterpret(LibXC.EnergyDensity{Cdouble}, expected[:, i])
                      for i in 1:size(expected, 2)],
                  [:ε, :v_a, :v_b, :δv_aa, :δv_ab, :δv_bb])
    elseif size(expected, 2) > 6
        cols = [:ε, :vrho_a, :vrho_b, :vsigma_aa, :vsigma_ab, :vsigma_bb, :v2rho_aa,
                :v2rho_ab, :v2rho_bb, :v2sigma2_aa_aa, :v2sigma2_aa_ab, :v2sigma2_aa_bb,
                :v2sigma2i_ab_ab, :v2sigma2_ab_bb, :v2sigma2_bb_bb, :v2rho_asigma_aa,
                :v2rho_asigma_ab, :v2rho_asigma_bb, :v2rho_bsigma_aa, :v2rho_bsigma_ab,
                :v2rho_bsigma_bb]
        result = DataFrame(Any[expected[:, i] for i in 1:size(expected, 2)], cols)
        cols = [:ε, :vrho_a, :vrho_b, :v2rho_aa, :v2rho_ab, :v2rho_bb]
        for u in cols
            result[u] = reinterpret(EnergyDensity{Cdouble}, result[u])
        end
    end
end

@testset "> LDA" begin
    input = input_data("BrOH")

    @testset ">> Unpolarizated " begin
        expected = expected_data("lda_x.BrOH.unpol.dat")
        ρ = input[:ρ_a] + input[:ρ_b]
        @test energy(:lda_x, ρ) ≈ expected[:ε]
        @test eltype(energy(:lda_x, ρ)) <: LibXC.EnergyDensity
        @test potential(:lda_x, ρ) ≈ expected[:v]
        @test eltype(potential(:lda_x, ρ)) <: LibXC.EnergyDensity
        @test second_energy_derivative(:lda_x, ρ) ≈ expected[:δv]
        @test eltype(second_energy_derivative(:lda_x, ρ)) <: LibXC.EnergyDensity

        func = XCFunctional(:lda_x, false)
        εxc, pot = energy_and_potential(func, ρ)
        @test εxc ≈ expected[:ε]
        @test pot ≈ expected[:v]

        εxc, pot, second_deriv, third_deriv = lda(:lda_x, ρ)
        @test εxc ≈ expected[:ε]
        @test pot ≈ expected[:v]
        @test second_deriv ≈ expected[:δv]
    end

    @testset ">> Polarized" begin
        expected = expected_data("lda_x.BrOH.pol.dat")

        ρs = vcat(input[:ρ_a]', input[:ρ_b]')
        @test energy(:lda_x, ρs) ≈ expected[:ε]
        @test potential(:lda_x, ρs) ≈ vcat(expected[:v_a]', expected[:v_b]')
        δv = vcat(expected[:δv_aa]', expected[:δv_ab]', expected[:δv_bb]')
        @test second_energy_derivative(:lda_x, ρs) ≈ δv

        func = XCFunctional(:lda_x, true)
        εxc, pot = energy_and_potential(func, ρs)
        @test εxc ≈ expected[:ε]
        @test pot ≈ vcat(expected[:v_a]', expected[:v_b]')

        εxc, pot, second_deriv, third_deriv = lda(:lda_x, ρs)
        @test εxc ≈ expected[:ε]
        @test pot ≈ vcat(expected[:v_a]', expected[:v_b]')
        @test second_deriv ≈ δv
    end
end

@testset "> GGA" begin
    input = input_data("BrOH")
#
#     @testset ">> Unpolarizated " begin
#         expected = expected_data("gga_c_pbe.BrOH.unpol.dat")
#         ρ = input[:ρ_a] + input[:ρ_b]
#         σ = input[:σ_aa] + 2input[:σ_ab] + input[:σ_bb]
#         @test_throws ArgumentError energy(:gga_c_pbe, ρ)
#         @test energy(:gga_c_pbe, ρ, σ) ≈ expected[:ε]
#
#         pot = potential(:gga_c_pbe, ρ, σ)
#         @test pot.rho ≈ expected[:v_a]
#         @test pot.sigma ≈ expected[:v_b]
#
#         second = second_energy_derivative(:gga_c_pbe, ρ, σ)
#         @test second.rho2 ≈ expected[:δv_aa]
#         @test second.rho_sigma ≈ expected[:δv_bb]
#         @test second.sigma2 ≈ expected[:δv_ab]
#
#         εxc, pot_rho, pot_sigma = energy_and_potential(:gga_c_pbe, ρ, σ)
#         @test εxc ≈ expected[:ε]
#         @test pot_rho ≈ expected[:v_a]
#         @test pot_sigma ≈ expected[:v_b]
#
#         all_out = gga(:gga_c_pbe, ρ, σ)
#         @test all_out.energy ≈ expected[:ε]
#         @test all_out.first_rho ≈ expected[:v_a]
#         @test all_out.first_sigma ≈ expected[:v_b]
#         @test all_out.second_rho2 ≈ expected[:δv_aa]
#         @test all_out.second_rho_sigma ≈ expected[:δv_bb]
#         @test all_out.second_sigma2 ≈ expected[:δv_ab]
#     end
#
#     @testset ">> Polarized" begin
#         expected = expected_data("gga_c_pbe.BrOH.pol.dat")
#
#         ρs = vcat(input[:ρ_a]', input[:ρ_b]')
#         σs = vcat(input[:σ_aa]', input[:σ_ab]', input[:σ_bb]')
#         @test energy(:gga_c_pbe, ρs, σs) ≈ expected[:ε]
#
#         pot = potential(:gga_c_pbe, ρs, σs)
#         @test pot.rho ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test pot.sigma ≈ expect
#
#         second = second_energy_derivative(:gga_c_pbe, ρs, σs)
#         v2rho = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         v2sigma = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
#                        expected[:v2sigma2_aa_bb]', expected[:v2sigma2i_ab_ab]',
#                        expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
#         v2rho = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         v2rho_sigma = vcat(expected[:v2rho_asigma_aa]', expected[:v2rho_asigma_ab]',
#                            expected[:v2rho_asigma_bb]', expected[:v2rho_bsigma_aa]',
#                            expected[:v2rho_bsigma_ab]', expected[:v2rho_bsigma_bb]')
#         @test second.rho2 ≈ v2rho
#         @test second.sigma2 ≈ v2sigma
#         @test second.rho_sigma ≈ v2rho_sigma
#
#         εxc, pot_rho, pot_sigma = energy_and_potential(:gga_c_pbe, ρs, σs)
#         @test εxc ≈ expected[:ε]
#         @test pot_rho ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test pot_sigma ≈ expect
#
#         all_out = gga(:gga_c_pbe, ρs, σs)
#         @test all_out.energy ≈ expected[:ε]
#         @test all_out.first_rho ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test all_out.first_sigma ≈ expect
#         v2rho = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         v2sigma = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
#                        expected[:v2sigma2_aa_bb]', expected[:v2sigma2i_ab_ab]',
#                        expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
#         v2rho = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         @test all_out.second_rho2 ≈ v2rho
#         @test all_out.second_rho_sigma ≈ v2rho_sigma
#         @test all_out.second_sigma2 ≈ v2sigma
#     end
end

end
