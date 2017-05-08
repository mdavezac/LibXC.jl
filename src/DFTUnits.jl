module DFTUnits
using Unitful
using UnitfulHartree

@unit ρ          "ρ"          Density                         1u"a₀^-3"      false
@unit ∇ρ         "∇ρ"         DensityGradient                 1u"a₀^-4"      false
@unit ∂ϵ_∂ρ      "∂ϵ_∂ρ"      FirstDensityDerivative          1u"Eₕ*a₀^3"    false
@unit ∂ϵ_∂∇ρ     "∂ϵ_∂∇ρ"     FirstGradientDerivative         1u"Eₕ*a₀^4"    false
@unit ∂²ϵ_∂ρ²    "∂²ϵ_∂ρ²"    SecondDensityDerivative         1u"Eₕ*a₀^6"  false
@unit ∂²ϵ_∂∇ρ²   "∂²ϵ_∂∇ρ²"   SecondGradientDerivative        1u"Eₕ*a₀^8"  false
@unit ∂²ϵ_∂ρ∂∇ρ  "∂²ϵ_∂ρ∂∇ρ"  SecondDensityGradientDerivative 1u"Eₕ*a₀^7"  false
@unit ∂³ϵ_∂ρ³    "∂³ϵ_∂ρ³"    ThirdDensityDerivative          1u"Eₕ*a₀^9"  false
@unit ∂³ϵ_∂∇ρ³   "∂³ϵ_∂∇ρ³"   ThirdGradientDerivative         1u"Eₕ*a₀^12" false
@unit ∂³ϵ_∂ρ²∂∇ρ "∂³ϵ_∂ρ²∂∇ρ" ThirdDensity2GradientDerivative 1u"Eₕ*a₀^10" false
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

@__dispatchers Density Ρ ρ "the electronic density"
@__dispatchers DensityGradient ∇Ρ ∇ρ "the gradient of the electronic density"
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
