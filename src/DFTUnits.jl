module DFTUnits
using Unitful
using UnitfulHartree

macro lintpragma(s) end

@lintpragma("Ignore use of undeclared variable Density")
@lintpragma("Ignore use of undeclared variable ρ")
@unit ρ          "ρ"          Density                         1u"a₀^-3"      false
@lintpragma("Ignore use of undeclared variable DensityGradient")
@lintpragma("Ignore use of undeclared variable ∇ρ")
@unit ∇ρ         "∇ρ"         DensityGradient                 1u"a₀^-4"      false
@lintpragma("Ignore use of undeclared variable FirstDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂ϵ_∂ρ")
@unit ∂ϵ_∂ρ      "∂ϵ_∂ρ"      FirstDensityDerivative          1u"Eₕ*a₀^3"    false
@lintpragma("Ignore use of undeclared variable FirstGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂ϵ_∂∇ρ")
@unit ∂ϵ_∂∇ρ     "∂ϵ_∂∇ρ"     FirstGradientDerivative         1u"Eₕ*a₀^4"    false
@lintpragma("Ignore use of undeclared variable SecondDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂ρ²")
@unit ∂²ϵ_∂ρ²    "∂²ϵ_∂ρ²"    SecondDensityDerivative         1u"Eₕ*a₀^6"  false
@lintpragma("Ignore use of undeclared variable SecondGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂∇ρ²")
@unit ∂²ϵ_∂∇ρ²   "∂²ϵ_∂∇ρ²"   SecondGradientDerivative        1u"Eₕ*a₀^8"  false
@lintpragma("Ignore use of undeclared variable SecondDensityGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂ρ∂∇ρ")
@unit ∂²ϵ_∂ρ∂∇ρ  "∂²ϵ_∂ρ∂∇ρ"  SecondDensityGradientDerivative 1u"Eₕ*a₀^7"  false
@lintpragma("Ignore use of undeclared variable ThirdDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ³")
@unit ∂³ϵ_∂ρ³    "∂³ϵ_∂ρ³"    ThirdDensityDerivative          1u"Eₕ*a₀^9"  false
@lintpragma("Ignore use of undeclared variable ThirdGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂∇ρ³")
@unit ∂³ϵ_∂∇ρ³   "∂³ϵ_∂∇ρ³"   ThirdGradientDerivative         1u"Eₕ*a₀^12" false
@lintpragma("Ignore use of undeclared variable ThirdDensity2GradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ²∂∇ρ")
@unit ∂³ϵ_∂ρ²∂∇ρ "∂³ϵ_∂ρ²∂∇ρ" ThirdDensity2GradientDerivative 1u"Eₕ*a₀^10" false
@lintpragma("Ignore use of undeclared variable ThirdDensityGradient2Derivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ∂∇ρ²")
@unit ∂³ϵ_∂ρ∂∇ρ² "∂³ϵ_∂ρ∂∇ρ²" ThirdDensityGradient2Derivative 1u"Eₕ*a₀^11" false

macro __dispatchers(name, abbr, units, docstr)
    dims = typeof(dimension(eval(units)))
    docstring = """ Dispatcher type over the dimension of $docstr """
    esc(quote
        @doc $docstring ->
        Unitful.Compat.@compat ($name){T, U} = Unitful.Quantity{T, $dims, U}
        @doc $docstring ->
        const $abbr = $name
    end)
end

@lintpragma("Ignore use of undeclared variable Ρ")
@__dispatchers Density Ρ ρ "the electronic density"
@lintpragma("Ignore use of undeclared variable ∇Ρ")
@__dispatchers DensityGradient ∇Ρ ∇ρ "the gradient of the electronic density"
@lintpragma("Ignore use of undeclared variable EnergyDensity")
@lintpragma("Ignore use of undeclared variable Ε")
@__dispatchers EnergyDensity Ε UnitfulHartree.Eₕ "the energy density"

const rho = ρ
const grho = ∇ρ
const ϵ = UnitfulHartree.Eₕ
const Exc = UnitfulHartree.Eₕ

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(DFTUnits)
end
end
