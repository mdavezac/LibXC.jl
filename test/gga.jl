using .Fixtures: input_data, expected_data, nospin, withspin
input = input_data("BrOH")

@testset ">> Unpolarizated " begin
    expected = expected_data("gga_c_pbe.BrOH.unpol.dat")
    expected[:v_b] = reinterpret(DHA.∂ϵ_∂σ{Cdouble}, ustrip(expected[:v_b]))
    expected[:δv_ab] = reinterpret(DHA.∂²ϵ_∂σ²{Cdouble}, ustrip(expected[:δv_ab]))
    expected[:δv_bb] = reinterpret(DHA.∂²ϵ_∂ρ∂σ{Cdouble}, ustrip(expected[:δv_bb]))

    ρ₀ = nospin(DHA.ρ{Cdouble}, input[:ρ_a] .+ input[:ρ_b])
    σ₀ = nospin(DHA.σ{Cdouble}, input[:σ_aa] .+ 2input[:σ_ab] .+ input[:σ_bb])

    functional = XCFunctional(:gga_c_pbe, false)

    @testset ">>> $name" for (name, (ρ, σ)) in [:AxisArray => (ρ₀, σ₀),
                                                :DenseArray => (ρ₀.data, σ₀.data)]
        @test_throws ArgumentError energy(functional, ρ)
        @test energy(functional, ρ, σ) ≈ expected[:ε]

        pot = potential(functional, ρ, σ)
        @test pot.∂ϵ_∂ρ ≈ expected[:v_a]
        @test pot.∂ϵ_∂σ ≈ expected[:v_b]

        second = second_energy_derivative(functional, ρ, σ)
        @test second.∂²ϵ_∂ρ² ≈ expected[:δv_aa]
        @test second.∂²ϵ_∂ρ∂σ ≈ expected[:δv_bb]
        @test second.∂²ϵ_∂σ² ≈ expected[:δv_ab]

        ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ = energy_and_potential(functional, ρ, σ)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ expected[:v_a]
        @test ∂ϵ_∂σ ≈ expected[:v_b]

        all_out = gga(functional, ρ, σ)
        @test all_out.ϵ ≈ expected[:ε]
        @test all_out.∂ϵ_∂ρ ≈ expected[:v_a]
        @test all_out.∂ϵ_∂σ ≈ expected[:v_b]
        @test all_out.∂²ϵ_∂ρ² ≈ expected[:δv_aa]
        @test all_out.∂²ϵ_∂ρ∂σ ≈ expected[:δv_bb]
        @test all_out.∂²ϵ_∂σ² ≈ expected[:δv_ab]
    end

    @testset ">>> Scalars" begin
        ρ, σ = 1.2u"ρ", 1.6u"σₑ"
        @test energy(functional, [ρ], [σ])[1] ≈ @inferred energy(:gga_c_pbe, ρ, σ)

        expected = potential(functional, [ρ], [σ])
        @test expected.∂ϵ_∂ρ[1] ≈ @inferred(potential(:gga_c_pbe, ρ, σ))[1]
        @test expected.∂ϵ_∂σ[1] ≈ @inferred(potential(:gga_c_pbe, ρ, σ))[2]

        expected = second_energy_derivative(functional, [ρ], [σ])
        @test expected.∂²ϵ_∂ρ²[1] ≈ @inferred(second_energy_derivative(:gga_c_pbe, ρ, σ))[1]
        @test expected.∂²ϵ_∂ρ∂σ[1] ≈ @inferred(second_energy_derivative(:gga_c_pbe, ρ, σ))[2]
        @test expected.∂²ϵ_∂σ²[1] ≈ @inferred(second_energy_derivative(:gga_c_pbe, ρ, σ))[3]

        expected = third_energy_derivative(XCFunctional(:gga_k_revapbe, false), [ρ], [σ])
        actual = @inferred third_energy_derivative(:gga_k_revapbe, ρ, σ)
        @test expected.∂³ϵ_∂ρ³[1] ≈ actual[1]
        @test expected.∂³ϵ_∂ρ²∂σ[1] ≈ actual[2]
        @test expected.∂³ϵ_∂ρ∂σ²[1] ≈ actual[3]
        @test expected.∂³ϵ_∂σ³[1] ≈ actual[4]
    end
end

@testset ">> Polarized" begin
    expected = expected_data("gga_c_pbe.BrOH.pol.dat")
    functional = XCFunctional(:gga_c_pbe, true)

    ρ₀ = withspin(DHA.ρ{Cdouble}, vcat(input[:ρ_a]', input[:ρ_b]'))
    σ₀ = withspin(DHA.σ{Cdouble}, vcat(input[:σ_aa]', input[:σ_ab]', input[:σ_bb]'))

    @testset ">>> $name" for (name, (ρs, σs)) in [:AxisArray => (ρ₀, σ₀),
                                                  :DenseArray => (ρ₀.data, σ₀.data)]
        @test energy(functional, ρs, σs) ≈ expected[:ε]
        rho = copy!(similar(ρs, typeof(1.0u"m^-3")), ρs)
        @test energy(functional, rho, σs) ≈ expected[:ε]

        pot = potential(functional, ρs, σs)
        @test pot.∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
        expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
        @test pot.∂ϵ_∂σ ≈ expect

        second = second_energy_derivative(functional, ρs, σs)
        ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
        ∂²ϵ_∂σ² = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
                       expected[:v2sigma2_aa_bb]', expected[:v2sigma2_ab_ab]',
                       expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
        ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
        ∂²ϵ_∂ρ∂σ = vcat(expected[:v2rho_asigma_aa]', expected[:v2rho_asigma_ab]',
                           expected[:v2rho_asigma_bb]', expected[:v2rho_bsigma_aa]',
                           expected[:v2rho_bsigma_ab]', expected[:v2rho_bsigma_bb]')
        @test second.∂²ϵ_∂ρ²  ≈ ∂²ϵ_∂ρ²
        @test second.∂²ϵ_∂σ²  ≈ ∂²ϵ_∂σ²
        @test second.∂²ϵ_∂ρ∂σ ≈ ∂²ϵ_∂ρ∂σ

        ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ = energy_and_potential(functional, ρs, σs)
        @test ϵ ≈ expected[:ε]
        @test ∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
        expect = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
        @test ∂ϵ_∂σ ≈ expect

        all_out = gga(functional, ρs, σs)
        @test all_out.ϵ ≈ expected[:ε]
        @test all_out.∂ϵ_∂ρ ≈ vcat(expected[:vrho_a]', expected[:vrho_b]')
        ∂ϵ_∂σ = vcat(expected[:vsigma_aa]', expected[:vsigma_ab]', expected[:vsigma_bb]')
        @test all_out.∂ϵ_∂σ ≈ ∂ϵ_∂σ
        ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
        ∂²ϵ_∂σ² = vcat(expected[:v2sigma2_aa_aa]', expected[:v2sigma2_aa_ab]',
                       expected[:v2sigma2_aa_bb]', expected[:v2sigma2_ab_ab]',
                       expected[:v2sigma2_ab_bb]', expected[:v2sigma2_bb_bb]')
        ∂²ϵ_∂ρ² = vcat(expected[:v2rho_aa]', expected[:v2rho_ab]', expected[:v2rho_bb]')
        @test all_out.∂²ϵ_∂ρ² ≈ ∂²ϵ_∂ρ²
        @test all_out.∂²ϵ_∂ρ∂σ ≈ ∂²ϵ_∂ρ∂σ
        @test all_out.∂²ϵ_∂σ² ≈ ∂²ϵ_∂σ²
    end

    @testset ">>> Scalars" begin
        ρα, ρβ, σαα, σαβ, σββ = 1.2u"ρ", 0.1u"ρ", 1.6u"σₑ", 1.2u"σₑ", 0.1u"σₑ" 
        expected = energy(functional, [ρα, ρβ], [σαα, σαβ, σββ])[1]
        @test expected ≈ @inferred energy(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)

        expected = potential(functional, [ρα, ρβ], [σαα, σαβ, σββ])
        @test expected.∂ϵ_∂ρ[1] ≈ @inferred(potential(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ))[1]
        @test expected.∂ϵ_∂ρ[2] ≈ potential(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)[2]
        @test expected.∂ϵ_∂σ[1] ≈ potential(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)[3]
        @test expected.∂ϵ_∂σ[2] ≈ potential(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)[4]
        @test expected.∂ϵ_∂σ[3] ≈ potential(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)[5]

        expected = second_energy_derivative(functional, [ρα, ρβ], [σαα, σαβ, σββ])
        actual = @inferred second_energy_derivative(:gga_c_pbe, ρα, ρβ, σαα, σαβ, σββ)
        @test expected.∂²ϵ_∂ρ² ≈ collect(actual[1:3])
        @test expected.∂²ϵ_∂ρ∂σ ≈ collect(actual[4:9])
        @test expected.∂²ϵ_∂σ² ≈ collect(actual[10:end])

        expected = third_energy_derivative(XCFunctional(:gga_k_revapbe, true),
                                           [ρα, ρβ], [σαα, σαβ, σββ])
        actual = third_energy_derivative(:gga_k_revapbe, ρα, ρβ, σαα, σαβ, σββ)
        @test expected.∂³ϵ_∂ρ³ ≈ collect(actual[1:4])
        @test expected.∂³ϵ_∂ρ²∂σ ≈ collect(actual[5:13])
        @test expected.∂³ϵ_∂ρ∂σ² ≈ collect(actual[14:25])
        @test expected.∂³ϵ_∂σ³ ≈ collect(actual[26:end])
    end
end

