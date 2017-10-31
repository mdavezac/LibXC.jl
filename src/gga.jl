module GGA

using AxisArrays
using StaticArrays
using ArgCheck

import ...LibXC: energy, energy!, first_energy_derivative, first_energy_derivative!,
                 energy_and_first_derivative, energy_and_first_derivative!,
                 second_energy_derivative, second_energy_derivative!,
                 third_energy_derivative, third_energy_derivative!, gga, gga!
using ..Internals
using ..Internals: AbstractLibXCFunctional, CFuncType, libxc
using ..OutputTuples
using ..Constants
using ..FunctionalMacros: @_all_wrapper_functionals, @_wrapper_functionals,
                          @_scalar_functional

using DFTShims: ColinearSpinFirst, Dispatch, is_spin_polarized, components, SpinDegenerate,
                SpinCategory
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

@_all_wrapper_functionals energy gga xc_gga_exc AxisArray ϵ
@_all_wrapper_functionals(energy_and_first_derivative, gga, xc_gga_exc_vxc,
                          GGAEnergyAndFirstDerivative, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
@_all_wrapper_functionals(first_energy_derivative, gga, xc_gga_vxc, GGAFirstDerivative,
                          ∂ϵ_∂ρ, ∂ϵ_∂σ)
@_all_wrapper_functionals(second_energy_derivative, gga, xc_gga_fxc, GGASecondDerivative,
                          ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ²)
@_all_wrapper_functionals(third_energy_derivative, gga, xc_gga_kxc, GGAThirdDerivative, ∂³ϵ_∂ρ³,
                          ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
@_wrapper_functionals(10, gga, xc_gga, GGATuple, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ,
                      ∂²ϵ_∂σ², ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
@_wrapper_functionals(10, gga, xc_gga, GGATuple, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ,
                      ∂²ϵ_∂σ²)
@_wrapper_functionals 10 gga xc_gga GGATuple ϵ ∂ϵ_∂ρ ∂ϵ_∂σ
@_wrapper_functionals 10 gga xc_gga GGATuple ϵ

gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DD.Arrays.ρ{Cdouble},
    σ::DD.Arrays.σ{Cdouble}) = begin
    family(func) ≠ Constants.gga && throw(ArgumentError("input function is not LDA"))
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
        gga!(func, ρ, σ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ∂σ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂σ²{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ³{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ²∂σ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂ρ∂σ²{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂³ϵ_∂σ³{Cdouble}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        gga!(func, ρ, σ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ∂σ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂σ²{Cdouble}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Cdouble}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Cdouble}, Spin))
    elseif Constants.exc ∈ f
        gga!(func, ρ, σ, ϵ)
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end

unsafe_energy(func::AbstractLibXCFunctional{Cdouble}, ρ::Cdouble, σ::Cdouble) = begin
    result = MVector{1, Cdouble}(0e0)
    ccall((:xc_gga_exc, libxc),
          Void, (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, ρ, σ, result)
    result[1]
end
unsafe_energy(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble, ρβ::Cdouble,
              σαα::Cdouble, σαβ::Cdouble, σββ::Cdouble) = begin
    result = MVector{1, Cdouble}(0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    cσ = SVector{3, Cdouble}(σαα, σαβ, σββ)
    ccall((:xc_gga_exc, libxc),
          Void, (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, cσ, result)
    result[1]
end
unsafe_first_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                               ρ::Cdouble, σ::Cdouble) = begin
    ∂ρ = MVector{1, Cdouble}(0e0)
    ∂σ = MVector{1, Cdouble}(0e0)
    ccall((:xc_gga_vxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, ρ, σ, ∂ρ, ∂σ)
    ∂ρ[1], ∂σ[1]
end
unsafe_first_energy_derivative(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble,
                               ρβ::Cdouble, σαα::Cdouble, σαβ::Cdouble,
                               σββ::Cdouble) = begin
    ∂ρ = MVector{2, Cdouble}(0e0, 0e0)
    ∂σ = MVector{3, Cdouble}(0e0, 0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    cσ = SVector{3, Cdouble}(σαα, σαβ, σββ)
    ccall((:xc_gga_vxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, cσ, ∂ρ, ∂σ)
    ∂ρ[1], ∂ρ[2], ∂σ[1], ∂σ[2], ∂σ[3]
end
unsafe_second_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                         ρ::Cdouble, σ::Cdouble) = begin
    ∂ρ = MVector{1, Cdouble}(0e0)
    ∂ρσ = MVector{1, Cdouble}(0e0)
    ∂σ = MVector{1, Cdouble}(0e0)
    ccall((:xc_gga_fxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, 1, ρ, σ, ∂ρ, ∂ρσ, ∂σ)
    ∂ρ[1], ∂ρσ[1], ∂σ[1]
end
unsafe_second_energy_derivative(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble,
                                ρβ::Cdouble, σαα::Cdouble, σαβ::Cdouble,
                                σββ::Cdouble) = begin
    ∂ρ = MVector{3, Cdouble}(0e0, 0e0, 0e0)
    ∂ρσ = MVector{6, Cdouble}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂σ = MVector{6, Cdouble}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    cσ = SVector{3, Cdouble}(σαα, σαβ, σββ)
    ccall((:xc_gga_fxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, 1, cρ, cσ, ∂ρ, ∂ρσ, ∂σ)
    (
        ∂ρ[1], ∂ρ[2], ∂ρ[3],
        ∂ρσ[1], ∂ρσ[2], ∂ρσ[3], ∂ρσ[4], ∂ρσ[5], ∂ρσ[6],
        ∂σ[1], ∂σ[2], ∂σ[3], ∂σ[4], ∂σ[5], ∂σ[6]
    )
end

unsafe_third_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                               ρ::Cdouble, σ::Cdouble) = begin
    ∂ρ³ = MVector{1, Cdouble}(0e0)
    ∂ρ²∂σ = MVector{1, Cdouble}(0e0)
    ∂ρ∂σ² = MVector{1, Cdouble}(0e0)
    ∂σ³ = MVector{1, Cdouble}(0e0)
    ccall((:xc_gga_kxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, ρ, σ, ∂ρ³, ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
    ∂ρ³[1], ∂ρ²∂σ[1], ∂ρ∂σ²[1], ∂σ³[1]
end

unsafe_third_energy_derivative(func::AbstractLibXCFunctional{Cdouble}, ρα::Cdouble,
                               ρβ::Cdouble, σαα::Cdouble, σαβ::Cdouble,
                               σββ::Cdouble) = begin
    ∂ρ³ = MVector{4, Cdouble}(0e0, 0e0, 0e0, 0e0)
    ∂ρ²∂σ = MVector{9, Cdouble}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂ρ∂σ² = MVector{12, Cdouble}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂σ³ = MVector{10, Cdouble}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Cdouble}(ρα, ρβ)
    cσ = SVector{3, Cdouble}(σαα, σαβ, σββ)
    ccall((:xc_gga_kxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, 1, cρ, cσ, ∂ρ³, ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
    ∂ρ³..., ∂ρ²∂σ..., ∂ρ∂σ²..., ∂σ³...
end

@_scalar_functional energy (ρ::ρ, σ::σ) (ϵ,) identity
@_scalar_functional energy (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ) (ϵ,) identity
@_scalar_functional first_energy_derivative (ρ::ρ, σ::σ) (∂ϵ_∂ρ, ∂ϵ_∂σ) tuple
@_scalar_functional(first_energy_derivative, (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ),
                    (∂ϵ_∂ρ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂ϵ_∂σ, ∂ϵ_∂σ), SGGAFirstDerivative)
@_scalar_functional second_energy_derivative (ρ::ρ, σ::σ) (∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ²) tuple
@_scalar_functional(second_energy_derivative, (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ),
                    (
                         ∂²ϵ_∂ρ², ∂²ϵ_∂ρ², ∂²ϵ_∂ρ²,
                         ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂ρ∂σ,
                         ∂²ϵ_∂σ², ∂²ϵ_∂σ², ∂²ϵ_∂σ², ∂²ϵ_∂σ², ∂²ϵ_∂σ², ∂²ϵ_∂σ²
                    ), GGASecondDerivatives)
@_scalar_functional(third_energy_derivative, (ρ::ρ, σ::σ),
                    (∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³), tuple)
@_scalar_functional(third_energy_derivative, (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ),
                    (
                        ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ³,
                        ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ,
                        ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ²∂σ,
                        ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ²,
                        ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂ρ∂σ²,
                        ∂³ϵ_∂σ³,  ∂³ϵ_∂σ³, ∂³ϵ_∂σ³, ∂³ϵ_∂σ³,  ∂³ϵ_∂σ³, ∂³ϵ_∂σ³,
                        ∂³ϵ_∂σ³,  ∂³ϵ_∂σ³, ∂³ϵ_∂σ³, ∂³ϵ_∂σ³
                    ), GGAThirdDerivatives)


end
