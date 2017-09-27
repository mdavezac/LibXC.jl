module GGA

using AxisArrays
using StaticArrays
using ArgCheck

import ...LibXC: energy, energy!, potential, potential!, energy_and_potential,
                 energy_and_potential!, second_energy_derivative, second_energy_derivative!,
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
@_all_wrapper_functionals(energy_and_potential, gga, xc_gga_exc_vxc, GGAEnergyAndPotential,
                          ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
@_all_wrapper_functionals potential gga xc_gga_vxc GGAPotential ∂ϵ_∂ρ ∂ϵ_∂σ
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

gga(func::AbstractLibXCFunctional{Float64}, ρ::DD.Arrays.ρ{Float64},
    σ::DD.Arrays.σ{Float64}) = begin
    family(func) ≠ Constants.gga && throw(ArgumentError("input function is not LDA"))
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
        gga!(func, ρ, σ, ϵ,
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
        gga!(func, ρ, σ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ²{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂ρ∂σ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂²ϵ_∂σ²{Float64}, Spin))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ, ϵ,
             similar(ρ, DH.Scalars.∂ϵ_∂ρ{Float64}, Spin),
             similar(ρ, DH.Scalars.∂ϵ_∂σ{Float64}, Spin))
    elseif Constants.exc ∈ f
        gga!(func, ρ, σ, ϵ)
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end

unsafe_energy(func::AbstractLibXCFunctional{Float64}, ρ::Float64, σ::Float64) = begin
    result = MVector{1, Float64}(0e0)
    ccall((:xc_gga_exc, libxc),
          Void, (Ptr{CFuncType}, Cint, Ref{Float64}, Ref{Float64}, Ptr{Float64}),
          func.c_ptr, 1, ρ, σ, result)
    result[1]
end
unsafe_energy(func::AbstractLibXCFunctional{Float64}, ρα::Float64, ρβ::Float64,
              σαα::Float64, σαβ::Float64, σββ::Float64) = begin
    result = MVector{1, Float64}(0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    cσ = SVector{3, Float64}(σαα, σαβ, σββ)
    ccall((:xc_gga_exc, libxc),
          Void, (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, cσ, result)
    result[1]
end
unsafe_potential(func::AbstractLibXCFunctional{Float64}, ρ::Float64, σ::Float64) = begin
    ∂ρ = MVector{1, Float64}(0e0)
    ∂σ = MVector{1, Float64}(0e0)
    ccall((:xc_gga_vxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, ρ, σ, ∂ρ, ∂σ)
    ∂ρ[1], ∂σ[1]
end
unsafe_potential(func::AbstractLibXCFunctional{Float64}, ρα::Float64, ρβ::Float64,
                 σαα::Float64, σαβ::Float64, σββ::Float64) = begin
    ∂ρ = MVector{2, Float64}(0e0, 0e0)
    ∂σ = MVector{3, Float64}(0e0, 0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    cσ = SVector{3, Float64}(σαα, σαβ, σββ)
    ccall((:xc_gga_vxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, cσ, ∂ρ, ∂σ)
    ∂ρ[1], ∂ρ[2], ∂σ[1], ∂σ[2], ∂σ[3]
end
unsafe_second_energy_derivative(func::AbstractLibXCFunctional{Float64},
                         ρ::Float64, σ::Float64) = begin
    ∂ρ = MVector{1, Float64}(0e0)
    ∂ρσ = MVector{1, Float64}(0e0)
    ∂σ = MVector{1, Float64}(0e0)
    ccall((:xc_gga_fxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}),
          func.c_ptr, 1, ρ, σ, ∂ρ, ∂ρσ, ∂σ)
    ∂ρ[1], ∂ρσ[1], ∂σ[1]
end
unsafe_second_energy_derivative(func::AbstractLibXCFunctional{Float64}, ρα::Float64,
                                ρβ::Float64, σαα::Float64, σαβ::Float64,
                                σββ::Float64) = begin
    ∂ρ = MVector{3, Float64}(0e0, 0e0, 0e0)
    ∂ρσ = MVector{6, Float64}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂σ = MVector{6, Float64}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    cσ = SVector{3, Float64}(σαα, σαβ, σββ)
    ccall((:xc_gga_fxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}),
          func.c_ptr, 1, cρ, cσ, ∂ρ, ∂ρσ, ∂σ)
    (
        ∂ρ[1], ∂ρ[2], ∂ρ[3],
        ∂ρσ[1], ∂ρσ[2], ∂ρσ[3], ∂ρσ[4], ∂ρσ[5], ∂ρσ[6],
        ∂σ[1], ∂σ[2], ∂σ[3], ∂σ[4], ∂σ[5], ∂σ[6]
    )
end

unsafe_third_energy_derivative(func::AbstractLibXCFunctional{Float64},
                               ρ::Float64, σ::Float64) = begin
    ∂ρ³ = MVector{1, Float64}(0e0)
    ∂ρ²∂σ = MVector{1, Float64}(0e0)
    ∂ρ∂σ² = MVector{1, Float64}(0e0)
    ∂σ³ = MVector{1, Float64}(0e0)
    ccall((:xc_gga_kxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, ρ, σ, ∂ρ³, ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
    ∂ρ³[1], ∂ρ²∂σ[1], ∂ρ∂σ²[1], ∂σ³[1]
end

unsafe_third_energy_derivative(func::AbstractLibXCFunctional{Float64}, ρα::Float64,
                               ρβ::Float64, σαα::Float64, σαβ::Float64,
                               σββ::Float64) = begin
    ∂ρ³ = MVector{4, Float64}(0e0, 0e0, 0e0, 0e0)
    ∂ρ²∂σ = MVector{9, Float64}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂ρ∂σ² = MVector{12, Float64}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    ∂σ³ = MVector{10, Float64}(0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0)
    cρ = SVector{2, Float64}(ρα, ρβ)
    cσ = SVector{3, Float64}(σαα, σαβ, σββ)
    ccall((:xc_gga_kxc, libxc),
          Void,
          (Ptr{CFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}),
          func.c_ptr, 1, cρ, cσ, ∂ρ³, ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
    ∂ρ³..., ∂ρ²∂σ..., ∂ρ∂σ²..., ∂σ³...
end

@_scalar_functional energy (ρ::ρ, σ::σ) (ϵ,) identity
@_scalar_functional energy (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ) (ϵ,) identity
@_scalar_functional potential (ρ::ρ, σ::σ) (∂ϵ_∂ρ, ∂ϵ_∂σ) tuple
@_scalar_functional(potential, (ρα::ρ, ρβ::ρ, σαα::σ, σβα::σ, σββ::σ),
                    (∂ϵ_∂ρ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂ϵ_∂σ, ∂ϵ_∂σ), GGAPotentials)
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
