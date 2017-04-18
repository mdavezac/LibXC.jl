module DFTUnits
using Unitful
using UnitfulHartree

Unitful.@dimension 𝐄 "𝐄" Electron # \mbfE
Unitful.@refunit 𝐞 "𝐞" Electronish 𝐄 false # \mbfe

@unit ρ         "ρ"         Density                           𝐞*1u"a₀^-3"        false
@unit ∇ρ         "∇ρ"         DensityGradient                 𝐞*1u"a₀^-4"        false
@unit ϵ         "ϵ"         EnergyDensity                     1u"Eₕ"/𝐞           false
@unit ∂ϵ_∂ρ     "∂ϵ_∂ρ"     FirstDensityDerivative            𝐞*1u"Eₕ^-2*a₀^3"   false
@unit ∂ϵ_∂∇ρ     "∂ϵ_∂∇ρ"     FirstGradientDerivative         𝐞*1u"Eₕ^-2*a₀^4"   false
@unit ∂²ϵ_∂ρ²   "∂²ϵ_∂ρ²"   SecondDensityDerivative           𝐞*1u"Eₕ^-3*a₀^6"   false
@unit ∂²ϵ_∂∇ρ²   "∂²ϵ_∂∇ρ²"   SecondGradientDerivative        𝐞*1u"Eₕ^-3*a₀^8"   false
@unit ∂²ϵ_∂ρ∂∇ρ  "∂²ϵ_∂ρ∂∇ρ"  SecondDensityGradientDerivative 𝐞*1u"Eₕ^-3*a₀^7"   false
@unit ∂³ϵ_∂ρ³   "∂³ϵ_∂ρ³"   ThirdDensityDerivative            𝐞*1u"Eₕ^-4*a₀^9"   false
@unit ∂³ϵ_∂∇ρ³   "∂³ϵ_∂∇ρ³"   ThirdGradientDerivative         𝐞*1u"Eₕ^-4*a₀^12"  false
@unit ∂³ϵ_∂ρ²∂∇ρ "∂³ϵ_∂ρ²∂∇ρ" ThirdDensity2GradientDerivative 𝐞*1u"Eₕ^-4*a₀^10"  false
@unit ∂³ϵ_∂ρ∂∇ρ² "∂³ϵ_∂ρ∂∇ρ²" ThirdDensityGradient2Derivative 𝐞*1u"Eₕ^-4*a₀^11"  false

const rho = ρ
const grho = ∇ρ
const Exc = ϵ

const localunits = Unitful.basefactors
const localpromotion = Unitful.promotion
function __init__()
    merge!(Unitful.basefactors, localunits)
    merge!(Unitful.promotion, localpromotion)
    Unitful.register(DFTUnits)
end
end
