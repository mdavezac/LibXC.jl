module GGA
using AxisArrays

import ...LibXC: energy, energy!, potential, potential!, energy_and_potential,
                 energy_and_potential!, second_energy_derivative, second_energy_derivative!,
                 third_energy_derivative, third_energy_derivative!, gga, gga!
using ..Internals
using ..OutputTuples
using ..Constants
using ..FunctionalMacros: @_all_wrapper_functionals, @_wrapper_functionals

using DFTShims: ColinearSpinFirst, Dispatch, is_spin_polarized, components, SpinDegenerate,
                SpinCategory
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

@_wrapper_functionals energy gga xc_gga_exc AxisArray ϵ
@_wrapper_functionals(energy_and_potential, gga, xc_gga_exc_vxc, GGAEnergyAndPotential,
                      ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
@_wrapper_functionals potential gga xc_gga_vxc GGAPotential ∂ϵ_∂ρ ∂ϵ_∂σ
@_wrapper_functionals(second_energy_derivative, gga, xc_gga_fxc, GGASecondDerivative,
                      ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ²)
@_wrapper_functionals(third_energy_derivative, gga, xc_gga_kxc, GGAThirdDerivative, ∂³ϵ_∂ρ³,
                      ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
@_wrapper_functionals(gga, gga, xc_gga, GGATuple, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ,
                      ∂²ϵ_∂σ², ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
@_all_wrapper_functionals(10, gga, xc_gga, GGATuple, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ,
                          ∂²ϵ_∂σ²)
@_all_wrapper_functionals 10 gga xc_gga GGATuple ϵ ∂ϵ_∂ρ ∂ϵ_∂σ
@_all_wrapper_functionals 10 gga xc_gga GGATuple ϵ

gga(func::AbstractLibXCFunctional{Float64}, ρ::DD.AxisArrays.ρ{Float64},
    σ::DD.AxisArrays.σ{Float64}) = begin
    family(func) ≠ Constants.gga && throw(ArgumentError("input function is not LDA"))
    Spin = SpinCategory(ρ)

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        gga!(func, ρ, σ,
             similar(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate()),
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂σ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ³{Float64}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ²∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ∂σ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂σ³{Float64}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        gga!(func, ρ, σ,
             similar(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate()),
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂σ²{Float64}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ,
             similar(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate()),
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Float64}, Spin))
    elseif Constants.exc ∈ f
        gga!(func, ρ, σ, similar(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate()))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
end
