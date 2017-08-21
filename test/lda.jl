using .Fixtures: input_data, expected_data, nospin, withspin

input = input_data("BrOH")

@testset ">> Unpolarizated " begin
    expected = expected_data("lda_x.BrOH.unpol.dat")
    ρ = nospin(DHA.ρ{Cdouble}, input[:ρ_a] + input[:ρ_b])
    @test energy(:lda_x, ρ) ≈ expected[:ε]
    @test eltype(energy(:lda_x, ρ)) <: DHA.ϵ
    @test potential(:lda_x, ρ) ≈ expected[:v]
    @test eltype(potential(:lda_x, ρ)) <: DHA.∂ϵ_∂ρ
    @test second_energy_derivative(:lda_x, ρ) ≈ expected[:δv]
    @test eltype(second_energy_derivative(:lda_x, ρ)) <: DHA.∂²ϵ_∂ρ²

    func = XCFunctional(:lda_x, false)
    ϵ, ∂ϵ_∂ρ = energy_and_potential(func, ρ)
    @test ϵ ≈ expected[:ε]
    @test ∂ϵ_∂ρ ≈ expected[:v]

    # checks unit and type conversion
    rho = uconvert(u"m^-3", ρ)
    with_conv = energy_and_potential(:lda_x, rho)
    @test with_conv[1] ≈ ϵ
    @test with_conv[2] ≈ ∂ϵ_∂ρ

    ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(:lda_x, ρ)
    @test ϵ ≈ expected[:ε]
    @test ∂ϵ_∂ρ ≈ expected[:v]
    @test ∂²ϵ_∂ρ² ≈ expected[:δv]
end

# @testset ">> Polarized" begin
#     expected = expected_data("lda_x.BrOH.pol.dat")
#
#     ρs = reinterpret(DHA.ρ{Cdouble}, vcat(input[:ρ_a]', input[:ρ_b]'))
#     @test energy(:lda_x, ρs) ≈ expected[:ε]
#     @test potential(:lda_x, ρs) ≈ vcat(expected[:v_a]', expected[:v_b]')
#     δv = vcat(expected[:δv_aa]', expected[:δv_ab]', expected[:δv_bb]')
#     @test second_energy_derivative(:lda_x, ρs) ≈ δv
#
#     func = XCFunctional(:lda_x, true)
#     ϵ, ∂ϵ_∂ρ = energy_and_potential(func, ρs)
#     @test ϵ ≈ expected[:ε]
#     @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')
#
#     ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(:lda_x, ρs)
#     @test ϵ ≈ expected[:ε]
#     @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')
#     @test ∂²ϵ_∂ρ² ≈ δv
# end
