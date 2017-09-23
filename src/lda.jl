module LDA

using AxisArrays
using StaticArrays
using ArgCheck
using Unitful

import ...LibXC: energy, energy!, potential, potential!, energy_and_potential,
                 energy_and_potential!, second_energy_derivative, second_energy_derivative!,
                 third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Internals
using ..Internals: AbstractLibXCFunctional, CFuncType, libxc
using ..OutputTuples
using ..Constants
using ..FunctionalMacros: @_all_wrapper_functionals, @_wrapper_functionals,
                          @_scalar_functional

using DFTShims: ColinearSpinFirst, Dispatch, is_spin_polarized, components, SpinDegenerate,
                SpinCategory
const DH = Dispatch.Hartree
const DD = Dispatch.Hartree


@_all_wrapper_functionals energy lda xc_lda_exc identity ϵ
@_all_wrapper_functionals energy_and_potential lda xc_lda_exc_vxc LDAEnergyAndPotential ϵ ∂ϵ_∂ρ
@_all_wrapper_functionals potential lda xc_lda_vxc identity ∂ϵ_∂ρ
@_all_wrapper_functionals second_energy_derivative lda xc_lda_fxc identity ∂²ϵ_∂ρ²
@_all_wrapper_functionals third_energy_derivative lda xc_lda_kxc identity ∂³ϵ_∂ρ³
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ² ∂³ϵ_∂ρ³
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ²
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ

lda(func::AbstractLibXCFunctional{Float64}, ρ::DD.Arrays.ρ{Float64}) = begin
    family(func) ≠ Constants.lda && throw(ArgumentError("input function is not LDA"))
    Spin = SpinCategory(func)

    # energy is a bit different since it is never spin polarized
    if ρ isa AxisArray || Spin isa SpinDegenerate
        ϵ = similar(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate())
    else
        # and axis less arrays can't deal with that easily
        @argcheck size(ρ, 1) == length(components(eltype(ρ), Spin))
        ϵ = similar(ρ, DH.Scalars.ϵ{Float64}, Base.tail(size(ρ)))
    end

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ³{Float64}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Float64}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin))
    elseif Constants.exc ∈ f
        lda!(func, ρ, ϵ)
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end

for (name, cname) in [:energy => :e, :potential => :v,
                      :second_energy_derivative => :f,
                      :third_energy_derivative => :k]
    @eval $name(func::AbstractLibXCFunctional{Float64}, ρ::Float64) = begin
        @argcheck spin(func) == Constants.unpolarized
        result = MVector{1, Float64}(0e0)
        ccall(($(QuoteNode(Symbol(:xc_lda_, cname, :xc))), libxc),
              Void, (Ptr{CFuncType}, Cint, Ref{Float64}, Ptr{Float64}),
              func.c_ptr, 1, ρ, result)
        result[1]
    end
end
energy(func::AbstractLibXCFunctional{Float64}, ρα::Float64, ρβ::Float64) = begin
    @argcheck spin(func) == Constants.polarized
    result = MVector{1, Float64}(0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    ccall((:xc_lda_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, result)
    result[1]
end
potential(func::AbstractLibXCFunctional{Float64}, ρα::Float64, ρβ::Float64) = begin
    @argcheck spin(func) == Constants.polarized
    result = MVector{2, Float64}(0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    ccall((:xc_lda_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, result)
    result
end
second_energy_derivative(func::AbstractLibXCFunctional{Float64},
                         ρα::Float64, ρβ::Float64) = begin
    @argcheck spin(func) == Constants.polarized
    result = MVector{3, Float64}(0e0, 0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    ccall((:xc_lda_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, result)
    result
end
third_energy_derivative(func::AbstractLibXCFunctional{Float64},
                         ρα::Float64, ρβ::Float64) = begin
    @argcheck spin(func) == Constants.polarized
    result = MVector{4, Float64}(0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    ccall((:xc_lda_kxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, result)
    result
end
energy_and_potential(func::AbstractLibXCFunctional{Float64}, ρ::Float64) = begin
    @argcheck spin(func) == Constants.unpolarized
    ϵ = MVector{1, Float64}(0e0)
    ∂ϵ_∂ρ = MVector{1, Float64}(0e0)
    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, ρ, ϵ, ∂ϵ_∂ρ)
    ϵ[1], ∂ϵ_∂ρ[1]
end
energy_and_potential(func::AbstractLibXCFunctional{Float64},
                     ρα::Float64, ρβ::Float64) = begin
    @argcheck spin(func) == Constants.polarized
    ϵ = MVector{1, Float64}(0e0)
    ∂ϵ_∂ρ = MVector{2, Float64}(0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, ϵ, ∂ϵ_∂ρ)
    ϵ[1], ∂ϵ_∂ρ[1], ∂ϵ_∂ρ[2]
end

@_scalar_functional energy (ρ::ρ,) (ϵ,)
@_scalar_functional energy (ρα::ρ, ρβ::ρ) (ϵ,)
@_scalar_functional potential (ρ::ρ,) (∂ϵ_∂ρ,)
@_scalar_functional potential (ρα::ρ, ρβ::ρ) (∂ϵ_∂ρ, ∂ϵ_∂ρ)
@_scalar_functional energy_and_potential (ρ::ρ,) (ϵ, ∂ϵ_∂ρ)
@_scalar_functional energy_and_potential (ρα::ρ, ρβ::ρ) (ϵ, ∂ϵ_∂ρ, ∂ϵ_∂ρ)
@_scalar_functional second_energy_derivative (ρ::ρ,) (∂²ϵ_∂ρ², )
@_scalar_functional second_energy_derivative (ρα::ρ, ρβ::ρ) (∂²ϵ_∂ρ², ∂²ϵ_∂ρ², ∂²ϵ_∂ρ²)
@_scalar_functional third_energy_derivative (ρ::ρ, ) (∂³ϵ_∂ρ³,)
@_scalar_functional(third_energy_derivative,
                    (ρα::ρ, ρβ::ρ), (∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³))

end
