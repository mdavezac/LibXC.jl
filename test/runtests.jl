module LibXCTests
using LibXC
using Base.Test
using Unitful
using DataFrames: DataFrame
using DFTShims

const DHA = Dispatch.Hartree.Scalars

include("fixtures.jl")
using .Fixtures

@testset "> Internal API" begin
    @test LibXC.LIB_VERSION == v"3.0.0"
    @test LibXC.FUNCTIONALS[:lda_x] == 1
    @test LibXC.FUNCTIONALS[:mgga_x_ms2] == 223
end

@testset "> Internal Functional API" begin
    include("internals.jl")
end

@testset "> Output tuples" begin
    include("OutputTuples.jl")
end

@testset "> LDA" begin
    include("lda.jl")
end

#
# @testset "> Array unit conversion"  begin
#     with_units = [1, 2, 4]u"m"
#     @test eltype(ustrip(LibXC.Units.conversion(u"cm", with_units))) == Cdouble
#     @test unit(eltype(LibXC.Units.conversion(u"cm", with_units))) == u"cm"
#     @test LibXC.Units.conversion(u"cm", with_units) ≈ with_units
#
#     @test LibXC.Units.conversion(u"cm", ustrip(with_units)) ≈ with_units // 100
#     @test eltype(ustrip(LibXC.Units.conversion(u"cm", ustrip(with_units)))) == Cdouble
#     @test unit(eltype(LibXC.Units.conversion(u"cm", ustrip(with_units)))) == u"cm"
#
#     @test unit(eltype(LibXC.Units.conversion(LibXC.Units.ϵ, [1, 2]))) == unit(LibXC.Units.ϵ)
#     @test eltype(ustrip(LibXC.Units.conversion(LibXC.Units.ϵ, [1, 2]))) == Cdouble
#     @test LibXC.Units.conversion(LibXC.Units.ϵ, [1, 2]) ≈ [1, 2] * LibXC.Units.ϵ{Cdouble}(1)
# end
#
# @testset "> LDA" begin
#     input = input_data("BrOH")
#
#     @testset ">> Unpolarizated " begin
#         expected = expected_data("lda_x.BrOH.unpol.dat")
#         ρ = reinterpret(LibXC.Units.ρ{Cdouble}, input[:ρ_a] + input[:ρ_b])
#         @test energy(:lda_x, ρ) ≈ expected[:ε]
#         @test eltype(energy(:lda_x, ρ)) <: LibXC.Units.ϵ
#         @test potential(:lda_x, ρ) ≈ expected[:v]
#         @test eltype(potential(:lda_x, ρ)) <: LibXC.Units.∂ϵ_∂ρ
#         @test second_energy_derivative(:lda_x, ρ) ≈ expected[:δv]
#         @test eltype(second_energy_derivative(:lda_x, ρ)) <: LibXC.Units.∂²ϵ_∂ρ²
#
#         func = XCFunctional(:lda_x, false)
#         ϵ, ∂ϵ_∂ρ = energy_and_potential(func, ρ)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ expected[:v]
#
#         # checks unit and type conversion
#         rho = reinterpret(Cdouble, LibXC.Units.conversion(u"m^-3", ρ))
#         rho = Float32[r for r in rho]u"m^-3"
#         with_conv = energy_and_potential(:lda_x, rho)
#         @test with_conv[1] ≈ ϵ
#         @test with_conv[2] ≈ ∂ϵ_∂ρ
#
#         ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(:lda_x, ρ)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ expected[:v]
#         @test ∂²ϵ_∂ρ² ≈ expected[:δv]
#     end
#
#     @testset ">> Polarized" begin
#         expected = expected_data("lda_x.BrOH.pol.dat")
#
#         ρs = reinterpret(LibXC.Units.ρ{Cdouble}, vcat(input[:ρ_a]', input[:ρ_b]'))
#         @test energy(:lda_x, ρs) ≈ expected[:ε]
#         @test potential(:lda_x, ρs) ≈ vcat(expected[:v_a]', expected[:v_b]')
#         δv = vcat(expected[:δv_aa]', expected[:δv_ab]', expected[:δv_bb]')
#         @test second_energy_derivative(:lda_x, ρs) ≈ δv
#
#         func = XCFunctional(:lda_x, true)
#         ϵ, ∂ϵ_∂ρ = energy_and_potential(func, ρs)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')
#
#         ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(:lda_x, ρs)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')
#         @test ∂²ϵ_∂ρ² ≈ δv
#     end
# end
#
# @testset "> GGA" begin
#     input = input_data("BrOH")
#
#     @testset ">> Unpolarizated " begin
#         expected = expected_data("gga_c_pbe.BrOH.unpol.dat")
#         expected[:v_b] = reinterpret(LibXC.Units.∂ϵ_∂∇ρ{Cdouble}, expected[:v_b])
#         expected[:δv_ab] = reinterpret(LibXC.Units.∂²ϵ_∂∇ρ²{Cdouble}, expected[:δv_ab])
#         expected[:δv_bb] = reinterpret(LibXC.Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}, expected[:δv_bb])
#         ρ = reinterpret(LibXC.Units.ρ{Cdouble}, input[:ρ_a] + input[:ρ_b])
#         ∇ρ = reinterpret(LibXC.Units.∇ρ{Cdouble}, input[:∇ρ_aa] + 2input[:∇ρ_ab] + input[:∇ρ_bb])
#         @test_throws ArgumentError energy(:gga_c_pbe, ρ)
#         @test energy(:gga_c_pbe, ρ, ∇ρ) ≈ expected[:ε]
#
#         pot = potential(:gga_c_pbe, ρ, ∇ρ)
#         @test pot.∂ϵ_∂ρ ≈ expected[:v_a]
#         @test pot.∂ϵ_∂∇ρ ≈ expected[:v_b]
#
#         second = second_energy_derivative(:gga_c_pbe, ρ, ∇ρ)
#         @test second.∂²ϵ_∂ρ² ≈ expected[:δv_aa]
#         @test second.∂²ϵ_∂ρ∂∇ρ ≈ expected[:δv_bb]
#         @test second.∂²ϵ_∂∇ρ² ≈ expected[:δv_ab]
#
#         ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ = energy_and_potential(:gga_c_pbe, ρ, ∇ρ)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ expected[:v_a]
#         @test ∂ϵ_∂∇ρ ≈ expected[:v_b]
#
#         all_out = gga(:gga_c_pbe, ρ, ∇ρ)
#         @test all_out.ϵ ≈ expected[:ε]
#         @test all_out.∂ϵ_∂ρ ≈ expected[:v_a]
#         @test all_out.∂ϵ_∂∇ρ ≈ expected[:v_b]
#         @test all_out.∂²ϵ_∂ρ² ≈ expected[:δv_aa]
#         @test all_out.∂²ϵ_∂ρ∂∇ρ ≈ expected[:δv_bb]
#         @test all_out.∂²ϵ_∂∇ρ² ≈ expected[:δv_ab]
#     end
#
#     @testset ">> Polarized" begin
#         expected = expected_data("gga_c_pbe.BrOH.pol.dat")
#
#         ρs = reinterpret(LibXC.Units.ρ{Cdouble}, vcat(input[:ρ_a]', input[:ρ_b]'))
#         ∇ρs = reinterpret(LibXC.Units.∇ρ{Cdouble},
#                          vcat(input[:∇ρ_aa]', input[:∇ρ_ab]', input[:∇ρ_bb]'))
#         @test energy(:gga_c_pbe, ρs, ∇ρs) ≈ expected[:ε]
#         @test energy(:gga_c_pbe, LibXC.Units.conversion(u"m^-3", ρs), ∇ρs) ≈ expected[:ε]
#
#         pot = potential(:gga_c_pbe, ρs, ∇ρs)
#         @test pot.∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test pot.∂ϵ_∂∇ρ ≈ expect
#
#         second = second_energy_derivative(:gga_c_pbe, ρs, ∇ρs)
#         ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         ∂²ϵ_∂∇ρ² = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
#                        expected[:v2sigma2_aa_bb]', expected[:v2sigma2_ab_ab]',
#                        expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
#         ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         ∂²ϵ_∂ρ∂∇ρ = vcat(expected[:v2rho_asigma_aa]', expected[:v2rho_asigma_ab]',
#                            expected[:v2rho_asigma_bb]', expected[:v2rho_bsigma_aa]',
#                            expected[:v2rho_bsigma_ab]', expected[:v2rho_bsigma_bb]')
#         @test second.∂²ϵ_∂ρ²  ≈ ∂²ϵ_∂ρ²
#         @test second.∂²ϵ_∂∇ρ²  ≈ ∂²ϵ_∂∇ρ²
#         @test second.∂²ϵ_∂ρ∂∇ρ ≈ ∂²ϵ_∂ρ∂∇ρ
#
#         ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ = energy_and_potential(:gga_c_pbe, ρs, ∇ρs)
#         @test ϵ ≈ expected[:ε]
#         @test ∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test ∂ϵ_∂∇ρ ≈ expect
#
#         all_out = gga(:gga_c_pbe, ρs, ∇ρs)
#         @test all_out.ϵ ≈ expected[:ε]
#         @test all_out.∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
#         ∂ϵ_∂∇ρ = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
#         @test all_out.∂ϵ_∂∇ρ ≈ ∂ϵ_∂∇ρ
#         ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         ∂²ϵ_∂∇ρ² = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
#                        expected[:v2sigma2_aa_bb]', expected[:v2sigma2_ab_ab]',
#                        expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
#         ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
#         @test all_out.∂²ϵ_∂ρ² ≈ ∂²ϵ_∂ρ²
#         @test all_out.∂²ϵ_∂ρ∂∇ρ ≈ ∂²ϵ_∂ρ∂∇ρ
#         @test all_out.∂²ϵ_∂∇ρ² ≈ ∂²ϵ_∂∇ρ²
#     end
# end

end
