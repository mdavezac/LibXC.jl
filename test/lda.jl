using .Fixtures: input_data, expected_data, nospin, withspin

input = input_data("BrOH")

@testset ">> Unpolarizated " begin
    expected = expected_data("lda_x.BrOH.unpol.dat")
    ρ₀ = nospin(DHA.ρ{Cdouble}, input[:ρ_a] + input[:ρ_b])
    functional = XCFunctional(:lda_x, false)

    @lintpragma("Ignore use of undeclared variable name")
    @testset ">>> $name" for (name, ρ) in [:AxisArray => ρ₀, :DenseArray => ρ₀.data]
        @test energy(functional, ρ) ≈ expected[:ε]
        @test eltype(energy(functional, ρ)) <: DHA.ϵ
        @test potential(functional, ρ) ≈ expected[:v]
        @test eltype(potential(functional, ρ)) <: DHA.∂ϵ_∂ρ
        @test second_energy_derivative(functional, ρ) ≈ expected[:δv]
        @test eltype(second_energy_derivative(functional, ρ)) <: DHA.∂²ϵ_∂ρ²

        ϵ, ∂ϵ_∂ρ = energy_and_potential(functional, ρ)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ expected[:v]

        ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(functional, ρ)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ expected[:v]
        @test ∂²ϵ_∂ρ² ≈ expected[:δv]

        # checks unit and type conversion
        rho = copy!(similar(ρ, typeof(1.0u"m^-3")), ρ)
        with_conv = energy_and_potential(functional, rho)
        @test with_conv[1] ≈ expected[:ε]
        @test with_conv[2] ≈ expected[:v]
    end
end

@testset ">> Polarized" begin
    expected = expected_data("lda_x.BrOH.pol.dat")
    ρ₀ˢ = withspin(DHA.ρ{Cdouble}, vcat(input[:ρ_a]', input[:ρ_b]'))
    functional = XCFunctional(:lda_x, true)

    @testset ">>> $name" for (name, ρs) in [:AxisArray => ρ₀ˢ, :DenseArray => ρ₀ˢ.data]
        @test energy(functional, ρs) ≈ expected[:ε]
        @test potential(functional, ρs) ≈ vcat(expected[:v_a]', expected[:v_b]')
        δv = vcat(expected[:δv_aa]', expected[:δv_ab]', expected[:δv_bb]')
        @test second_energy_derivative(functional, ρs) ≈ δv

        ϵ, ∂ϵ_∂ρ = energy_and_potential(functional, ρs)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')

        ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³ = lda(functional, ρs)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ vcat(expected[:v_a]', expected[:v_b]')
        @test ∂²ϵ_∂ρ² ≈ δv
    end
end
