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

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable energy")
@lintpragma("Ignore use of undeclared variable xc_lda_exc")
@lintpragma("Ignore use of undeclared variable ϵ")
@lintpragma("Ignore use of undeclared variable SLDAEnergyAndPotential")
@lintpragma("Ignore use of undeclared variable energy_and_potential")
@lintpragma("Ignore use of undeclared variable xc_lda_exc_vxc")
@lintpragma("Ignore use of undeclared variable ∂ϵ_∂ρ")
@lintpragma("Ignore use of undeclared variable potential")
@lintpragma("Ignore use of undeclared variable xc_lda_vxc")
@lintpragma("Ignore use of undeclared variable xc_lda_fxc")
@lintpragma("Ignore use of undeclared variable xc_lda_kxc")
@lintpragma("Ignore use of undeclared variable second_energy_derivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂ρ²")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ³")
@lintpragma("Ignore use of undeclared variable third_energy_derivative")
@lintpragma("Ignore use of undeclared variable LDATuple")
@lintpragma("Ignore use of undeclared variable LDAEnergyAndPotential")
@lintpragma("Ignore use of undeclared variable SLDASecondDerivative")
@lintpragma("Ignore use of undeclared variable SLDAThirdDerivative")
@lintpragma("Ignore use of undeclared variable xc_lda")
@lintpragma("Ignore use of undeclared variable family")
@lintpragma("Ignore use of undeclared variable flags")
@lintpragma("Ignore use of undeclared variable lda!")
@lintpragma("Ignore use of undeclared variable ccall")
@lintpragma("Ignore use of undeclared variable ρ")
@lintpragma("Ignore use of undeclared variable ρα")
@lintpragma("Ignore use of undeclared variable ρβ")


@_all_wrapper_functionals energy lda xc_lda_exc identity ϵ
@_all_wrapper_functionals(energy_and_potential, lda, xc_lda_exc_vxc, LDAEnergyAndPotential,
                          ϵ, ∂ϵ_∂ρ)
@_all_wrapper_functionals potential lda xc_lda_vxc identity ∂ϵ_∂ρ
@_all_wrapper_functionals second_energy_derivative lda xc_lda_fxc identity ∂²ϵ_∂ρ²
@_all_wrapper_functionals third_energy_derivative lda xc_lda_kxc identity ∂³ϵ_∂ρ³
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ² ∂³ϵ_∂ρ³
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ ∂²ϵ_∂ρ²
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ ∂ϵ_∂ρ
@_wrapper_functionals 4 lda xc_lda LDATuple ϵ

lda(func::AbstractLibXCFunctional{Cdouble}, ρ::DD.Arrays.ρ{Cdouble}) = begin
    family(func) ≠ Constants.lda && throw(ArgumentError("input function is not LDA"))
    Spin = SpinCategory(func)

    # energy is a bit different since it is never spin polarized
    if ρ isa AxisArray || Spin isa SpinDegenerate
        ϵ = similar(ρ, DH.Scalars.ϵ{Cdouble}, SpinDegenerate())
    else
        # and axis less arrays can't deal with that easily
        @argcheck size(ρ, 1) == length(components(eltype(ρ), Spin))
        ϵ = similar(ρ, DH.Scalars.ϵ{Cdouble}, Base.tail(size(ρ)))
    end

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ³{Cdouble}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Cdouble}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        lda!(func, ρ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin))
    elseif Constants.exc ∈ f
        lda!(func, ρ, ϵ)
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end

for (name, cname) in [:unsafe_energy => :e, :unsafe_potential => :v,
                      :unsafe_second_energy_derivative => :f,
                      :unsafe_third_energy_derivative => :k]
    @eval $name(func::AbstractLibXCFunctional{Cdouble}, ρ::Cdouble) = begin
        result = MVector{1, Cdouble}(0e0)
        ccall(($(QuoteNode(Symbol(:xc_lda_, cname, :xc))), libxc),
              Void, (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ptr{Cdouble}),
              func.c_ptr, 1, ρ, result)
        result[1]
    end
end
unsafe_energy(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble, ρβ::Cdouble) = begin
    result = MVector{1, Cdouble}(0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    ccall((:xc_lda_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, result)
    result[1]
end
unsafe_potential(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble, ρβ::Cdouble) = begin
    result = MVector{2, Cdouble}(0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    ccall((:xc_lda_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, result)
    result
end
unsafe_second_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                                ρα::Cdouble, ρβ::Cdouble) = begin
    result = MVector{3, Cdouble}(0e0, 0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    ccall((:xc_lda_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, result)
    result
end
unsafe_third_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                               ρα::Cdouble, ρβ::Cdouble) = begin
    result = MVector{4, Cdouble}(0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    ccall((:xc_lda_kxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, result)
    result
end
unsafe_energy_and_potential(func::AbstractLibXCFunctional{Cdouble}, ρ::Cdouble) = begin
    ϵ = MVector{1, Cdouble}(0e0)
    ∂ϵ_∂ρ = MVector{1, Cdouble}(0e0)
    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, ρ, ϵ, ∂ϵ_∂ρ)
    ϵ[1], ∂ϵ_∂ρ[1]
end
unsafe_energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                            ρα::Cdouble, ρβ::Cdouble) = begin
    ϵ = MVector{1, Cdouble}(0e0)
    ∂ϵ_∂ρ = MVector{2, Cdouble}(0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, ϵ, ∂ϵ_∂ρ)
    ϵ[1], ∂ϵ_∂ρ[1], ∂ϵ_∂ρ[2]
end

@_scalar_functional energy (ρ::ρ,) (ϵ,) identity
@_scalar_functional energy (ρα::ρ, ρβ::ρ) (ϵ,) identity
@_scalar_functional potential (ρ::ρ,) (∂ϵ_∂ρ,) identity
@_scalar_functional potential (ρα::ρ, ρβ::ρ) (∂ϵ_∂ρ, ∂ϵ_∂ρ) tuple
@_scalar_functional energy_and_potential (ρ::ρ,) (ϵ, ∂ϵ_∂ρ) LDAEnergyAndPotential
@_scalar_functional(energy_and_potential, (ρα::ρ, ρβ::ρ), (ϵ, ∂ϵ_∂ρ, ∂ϵ_∂ρ),
                    SLDAEnergyAndPotential)
@_scalar_functional second_energy_derivative (ρ::ρ,) (∂²ϵ_∂ρ², ) identity
@_scalar_functional(second_energy_derivative, (ρα::ρ, ρβ::ρ),
                    (∂²ϵ_∂ρ², ∂²ϵ_∂ρ², ∂²ϵ_∂ρ²), SLDASecondDerivative)
@_scalar_functional third_energy_derivative (ρ::ρ, ) (∂³ϵ_∂ρ³,) identity
@_scalar_functional(third_energy_derivative,
                    (ρα::ρ, ρβ::ρ), (∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³),
                    SLDAThirdDerivative)

end
