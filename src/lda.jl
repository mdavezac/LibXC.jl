module LDA

using AxisArrays
import ...LibXC: energy, energy!, potential, potential!, energy_and_potential,
                 energy_and_potential!, second_energy_derivative, second_energy_derivative!,
                 third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Internals
using ..OutputTuples
using ..Constants
using ..FunctionalMacros: @_all_wrapper_functionals, @_wrapper_functionals

using DFTShims: ColinearSpinFirst, Dispatch, is_spin_polarized, components
const DH = Dispatch.Hartree
const DD = Dispatch.Hartree


@_wrapper_functionals energy lda xc_lda_exc AxisArray ϵ
@_wrapper_functionals energy_and_potential lda xc_lda_exc_vxc LDAEnergyAndPotential ϵ ∂ϵ_∂ρ
@_wrapper_functionals potential lda xc_lda_vxc AxisArray ∂ϵ_∂ρ
@_wrapper_functionals second_energy_derivative lda xc_lda_fxc AxisArray ∂²ϵ_∂ρ²
@_wrapper_functionals third_energy_derivative lda xc_lda_kxc AxisArray ∂³ϵ_∂ρ³
@_wrapper_functionals lda lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ² ∂³ϵ_∂ρ³
@_all_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ² 
@_all_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ
@_all_wrapper_functionals 4 lda xc_lda LDATuple ϵ

lda(func::AbstractLibXCFunctional{Float64}, ρ::DD.AxisArrays.ρ{Float64}) = begin
    family(func) ≠ Constants.lda && throw(ArgumentError("input function is not LDA"))

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        lda!(func, ρ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ²{Float64}, ρ),
             similar(DH.Scalars.∂³ϵ_∂ρ³{Float64}, ρ))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        lda!(func, ρ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ²{Float64}, ρ))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        lda!(func, ρ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ))
    elseif Constants.exc ∈ f
        lda!(func, ρ, similar(DH.Scalars.ϵ{Float64}, ρ))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end

end
