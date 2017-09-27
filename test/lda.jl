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

    @testset ">>> Scalars" begin
        ρ = ρ₀
        rho = copy!(similar(ρ, typeof(Float32(1.0)u"m^-3")), ρ)

        @test energy.(:lda_x, ρ) ≈ expected[:ε]
        @test energy.(:lda_x, rho) ≈ expected[:ε]
        @inferred energy(:lda_x, 1u"ρ")
        @inferred energy(:lda_x, 1u"nm^-3")

        @test potential.(:lda_x, ρ) ≈ expected[:v]
        @test potential.(:lda_x, rho) ≈ expected[:v]

        ρ = 1.0u"ρ"
        expected = energy_and_potential(functional, [ρ])
        @test expected.ϵ[1] == @inferred(energy_and_potential(:lda_x, ρ))[1]
        @test expected.∂ϵ_∂ρ[1] == @inferred(energy_and_potential(:lda_x, ρ))[2]

        expected = second_energy_derivative(functional, [ρ])[1]
        @test expected == @inferred second_energy_derivative(:lda_x, ρ)

        expected = third_energy_derivative(functional, [ρ])[1]
        @test expected == @inferred third_energy_derivative(:lda_x, ρ)
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

    @testset ">>> Scalars" begin
        ρα = view(ρ₀ˢ, Axis{:spin}(1))
        ρβ = view(ρ₀ˢ, Axis{:spin}(2))
        rhoα = copy!(similar(ρα, typeof(Float32(1.0)u"m^-3")), ρα)
        rhoβ = copy!(similar(ρβ, typeof(Float32(1.0)u"m^-3")), ρβ)

        @test energy.(functional, ρα, ρβ) ≈ expected[:ε]
        @test energy.(functional, rhoα, rhoβ) ≈ expected[:ε]
        @inferred energy(functional, 1u"ρ", 1u"nm^-3")
        @inferred energy(functional, 1u"nm^-3", 1.0u"ρ")

        a, b = 1.0u"ρ", 1.5u"ρ"
        expected = tuple(potential(XCFunctional(:lda_x, true), [a, b])...)
        @test expected == @inferred potential(:lda_x, a, b)

        expected = energy_and_potential(XCFunctional(:lda_x, true), [a, b])
        @test expected.ϵ[1] == @inferred(energy_and_potential(:lda_x, a, b))[1]
        @test expected.∂ϵ_∂ρ[1] == @inferred(energy_and_potential(:lda_x, a, b))[2]
        @test expected.∂ϵ_∂ρ[2] == @inferred(energy_and_potential(:lda_x, a, b))[3]

        expected = SLDASecondDerivative(second_energy_derivative(functional, [a, b])...)
        @test expected == @inferred second_energy_derivative(:lda_x, a, b)

        expected = SLDAThirdDerivative(third_energy_derivative(functional, [a, b])...)
        @test expected == @inferred third_energy_derivative(:lda_x, a, b)
    end
end
