module GGA
import LibXC: energy, energy!, potential, potential!, energy_and_potential,
              energy_and_potential!, second_energy_derivative, second_energy_derivative!,
              third_energy_derivative, third_energy_derivative!, gga, gga!
using ..Constants
using ..Internals
using ..Internals: CFuncType, libxc, XCFunctional
using ..OutputTuples
using ..Checks
using ..LDA: _ddargument, _dhargument, _deargument, _conversion, _convert_back,
             _convert_to_array, _check_availability, _check_axes_nospin, _check_axes_spin,
             _similar

using AxisArrays
using DocStringExtensions

using DFTShims: SpinCategory, SpinDegenerate, ColinearSpinFirst, Dispatch, is_spin_polarized
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable ccall")

const TYPES = (:ϵ, :∂ϵ_∂ρ, :∂ϵ_∂σ, :∂²ϵ_∂ρ², :∂²ϵ_∂ρ∂σ, :∂²ϵ_∂σ², :∂³ϵ_∂ρ³, :∂³ϵ_∂ρ²∂σ,
               :∂³ϵ_∂ρ∂σ², :∂³ϵ_∂σ³)
const FUNCTIONS = Dict(:energy => TYPES[1:1], :energy_and_potential => TYPES[1:3],
                       :potential => TYPES[2:3], :second_energy_derivative => TYPES[4:6],
                       :third_energy_derivative => TYPES[6:end], :gga => TYPES)
const OUTPUTS = Dict(:energy => :AxisArray, :energy_and_potential => :GGAEnergyAndPotential,
                     :potential => :GGAPotential,
                     :second_energy_derivative => :GGASecondDerivative,
                     :third_energy_derivative => :GGAThirdDerivative, :gga => :GGATuple)
const FUNCNAMES = Dict(:energy => :xc_gga_exc, :energy_and_potential => :xc_gga_exc_vxc,
                       :potential => :xc_gga_vxc, :second_energy_derivative => :xc_gga_fxc,
                       :third_energy_derivative => :xc_gga_kxc, :gga => :xc_gga)
const FUNCTYPES = Dict(:energy => :exc, :energy_and_potential => (:exc, :vxc),
                       :potential => :vxc, :second_energy_derivative => :fxc,
                       :third_energy_derivative => :kxc, :gga => (:exc, :vxc, :fxc, :kxc))

_check_functional(functype::Symbol) = begin
    msg = "Functional is GGA, but input array do not correspond"
    :(family($functype) == Constants.gga || throw(ArgumentError($msg)))
end
_empty_arrays(name::Symbol) = :(zeros(DH.Scalars.$name{Float64}, SpinCategory(ρ), (0, )))
_result_slice(name::Symbol) = :(result.$name)

for (func, outputs) in FUNCTIONS
    allargs = (:ρ, :σ, outputs...)
    otype = OUTPUTS[func]
    func! = Symbol(func, :!)
    @eval begin
        $func!(func::AbstractLibXCFunctional, $(_ddargument.(allargs)...)) = begin
            result = $func!(SpinCategory(func), func, $(_conversion.(allargs)...))
            $otype($(_convert_back(outputs)...))
        end
        $func!(::SpinDegenerate, func::AbstractLibXCFunctional,
               $(_dhargument.(allargs)...)) = begin
            $(_check_axes_nospin.(outputs)...)
            $(_check_functional(:func))
            $(_check_availability.(FUNCTYPES[func], :func))
            
            $func!(func, $(_convert_to_array.(allargs)...))

            $otype($(outputs...))
        end
        $func!(::ColinearSpinFirst, func::AbstractLibXCFunctional,
               $(_dhargument.(allargs)...)) = begin
            $(_check_axes_spin.(allargs)...)
            $(_check_functional(:func))
            $(_check_availability.(:func, FUNCTYPES[func]))

            $func!(func, $(_convert_to_array.(allargs)...))

            $otype($(outputs...))
        end

        $func!(func::AbstractLibXCFunctional, $(_deargument.(allargs)...)) = begin
            @lintpragma("Ignore use of undeclared variable ρ")
            @lintpragma("Ignore unused u")
            ccall(($(QuoteNode(FUNCNAMES[func])), $libxc), Void,
                  (Ptr{CFuncType}, Cint, $([:(Ptr{Float64}) for u in allargs]...)),
                  func.c_ptr, length(ρ) / (spin(func) == Constants.polarized ? 2: 1),
                  $(allargs...))
            $outputs
        end

        $func(func::AbstractLibXCFunctional, ρ::DD.AxisArrays.ρ, σ::DD.AxisArrays.σ) =
            $func!(func, ρ, σ, $(_similar.(outputs)...))

        $func(name::Symbol, ρ::DD.AxisArrays.ρ, σ::DD.AxisArrays.σ) =
            $func(XCFunctional(name, is_spin_polarized(ρ)), ρ, σ)
    end
end

# for i in (1, 3, 6)
#     allargs = (:ρ, :σ, TYPES[1:i]...)
#     @eval gga!(func::AbstractLibXCFunctional, $(_ddargument.(allargs)...)) = begin
#         result = gga!(SpinCategory(func), func, $(_conversion.(allargs)...))
#         @lintpragma("Ignore use of undeclared variable GGATuple")
#         GGATuple($(_convert_back(outputs)...))
#         # GGATuple($(_convert_back(TYPES[1:i)...)...)
#     end
# end

gga(func::AbstractLibXCFunctional{Float64}, ρ::DD.AxisArrays.ρ{Float64},
    σ::DD.AxisArrays.σ{Float64}) = begin
    @check_functional func gga

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        gga!(func, ρ, σ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂σ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ²{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ∂σ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂σ²{Float64}, ρ),
             similar(DH.Scalars.∂³ϵ_∂ρ³{Float64}, ρ),
             similar(DH.Scalars.∂³ϵ_∂ρ²∂σ{Float64}, ρ),
             similar(DH.Scalars.∂³ϵ_∂ρ∂σ²{Float64}, ρ),
             similar(DH.Scalars.∂³ϵ_∂σ³{Float64}, ρ))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        gga!(func, ρ, σ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂σ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ²{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂ρ∂σ{Float64}, ρ),
             similar(DH.Scalars.∂²ϵ_∂σ²{Float64}, ρ))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ,
             similar(DH.Scalars.ϵ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂ρ{Float64}, ρ),
             similar(DH.Scalars.∂ϵ_∂σ{Float64}, ρ))
    elseif Constants.exc ∈ f
        gga!(func, ρ, σ, similar(DH.Scalars.ϵ{Float64}, ρ))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
end
# Mostly so linting doesn't throw false positives when processing this file alone
# using LibXC: Constants, CFuncType, libxc, spin, output_size, flags
# using LibXC: GGAPotential, GGASecondDerivative, GGAThirdDerivative
# using LibXC: GGAEnergyAndPotential, AllGGA
# using LibXC.Checks: @check_functional, @check_availability, @check_size
# @lintpragma("Ignore use of undeclared variable ccall")
#
# """
#     $(SIGNATURES)
#
# GGA energy as a function of the density ρ and the contracted gradient density σ.
# When ρ is spin-polarized, σ=(∇ρ↑⋅∇ρ↑, ∇ρ↑⋅∇ρ↓, ∇ρ↓⋅∇ρ↓). The dimensionality is as follows:
#
# |GGA | unpolarized | polarized                |
# |----|-------------|--------------------------|
# |ρ   | any         | `(2, ...)`               |
# |σ   | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
# |ϵ   | `size(ρ)`   | `size(ρ)[2:end]`         |
#
# In other words, in spin-polarized systems, the components arising from the spin are the
# fastest changing direction (Julia is column-major).
# """
# function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}},
#                  σ::DenseArray{Units.σ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}})
#     energy!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ), reinterpret(Cdouble, ϵ))
#     ϵ
# end
# function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
#                  σ::DenseArray{Cdouble}, ϵ::DenseArray{Cdouble})
#     @check_functional func gga
#     @check_availability func exc
#     @check_size func ρ σ 3
#     @check_size func ρ ϵ 1
#
#     ccall((:xc_gga_exc, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
#           func.c_ptr, length(ρ) /convert(Int64, spin(func)), ρ, σ, ϵ)
#     ϵ
# end
# function energy{Ρ <: Quantity, ∇Σ <: Quantity}(func::AbstractLibXCFunctional,
#                                                ρ::DenseArray{Ρ}, σ::DenseArray{∇Σ})
#     energy!(func, Units.conversion(Units.ρ, ρ), Units.conversion(Units.σ, σ),
#             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
# end
# function energy(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
#                 σ::DenseArray{Cdouble})
#     energy!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
# end
#
#
# """
#     $(SIGNATURES)
#
# GGA potential computed in place. The dimensionality of the different arrays are as follows:
#
# |GGA     | unpolarized | polarized                |
# |--------|-------------|--------------------------|
# |ρ       | any         | `(2, ...)`               |
# |σ      | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
# |∂ϵ/∂ρ   | `size(ρ)`   | `size(ρ)`                |
# |∂ϵ/∂σ  | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
# """
# function potential!(func::AbstractLibXCFunctional{Cdouble},
#                     ρ::DenseArray{Units.ρ{Cdouble}},
#                     σ::DenseArray{Units.σ{Cdouble}},
#                     ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
#                     ∂ϵ_∂σ::DenseArray{Units.∂ϵ_∂σ{Cdouble}})
#     potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                reinterpret(Cdouble, ∂ϵ_∂ρ), reinterpret(Cdouble, ∂ϵ_∂σ))
#     GGAPotential(∂ϵ_∂ρ, ∂ϵ_∂σ)
# end
# function potential!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
#                     σ::DenseArray{Cdouble}, ∂ϵ_∂ρ::DenseArray{Cdouble},
#                     ∂ϵ_∂σ::DenseArray{Cdouble})
#     @check_functional func gga
#     @check_availability func vxc
#     @check_size func ρ σ 3
#     @check_size func ρ ∂ϵ_∂ρ 2
#     @check_size func ρ ∂ϵ_∂σ 3
#
#     ccall((:xc_gga_vxc, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
#           func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
#     GGAPotential(∂ϵ_∂ρ, ∂ϵ_∂σ)
# end
# function potential{Ρ <: Quantity, ∇Σ <: Quantity}(func::AbstractLibXCFunctional,
#                                                   ρ::DenseArray{Ρ},
#                                                   σ::DenseArray{∇Σ})
#     potential!(func, Units.conversion(Units.ρ, ρ), Units.conversion(Units.σ, σ),
#                similar(ρ, Units.∂ϵ_∂ρ{Cdouble}),
#                similar(ρ, Units.∂ϵ_∂σ{Cdouble}, output_size(func, ρ, 3)))
# end
# function potential(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
#                    σ::DenseArray{Cdouble})
#     potential!(func, ρ, σ, similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
# end
#
# """
#     $(SIGNATURES)
#
# Second derivatives of GGA energy w.r.t. ρ and σ. The dimensionality of the arrays is as
# follows:
#
# |GGA       | unpolarized | polarized                |
# |----------|-------------|--------------------------|
# |ρ         | any         | `(2, ...)`               |
# |σ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
# |∂²ϵ/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
# |∂²ϵ/∂ρ∂σ  | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
# |∂²ϵ/∂σ²   | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
# """
# function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
#                                    ρ::DenseArray{Units.ρ{Cdouble}},
#                                    σ::DenseArray{Units.σ{Cdouble}},
#                                    ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
#                                    ∂²ϵ_∂ρ∂σ::DenseArray{Units.∂²ϵ_∂ρ∂σ{Cdouble}},
#                                    ∂²ϵ_∂σ²::DenseArray{Units.∂²ϵ_∂σ²{Cdouble}})
#     second_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                               reinterpret(Cdouble, ∂²ϵ_∂ρ²), reinterpret(Cdouble, ∂²ϵ_∂ρ∂σ),
#                               reinterpret(Cdouble, ∂²ϵ_∂σ²))
#     GGASecondDerivative(∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ²)
# end
# function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
#                                    ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
#                                    ∂²ϵ_∂ρ²::DenseArray{Cdouble},
#                                    ∂²ϵ_∂ρ∂σ::DenseArray{Cdouble},
#                                    ∂²ϵ_∂σ²::DenseArray{Cdouble})
#     @check_functional func gga
#     @check_availability func fxc
#     @check_size func ρ σ 3
#     @check_size func ρ ∂²ϵ_∂ρ² 3
#     @check_size func ρ ∂²ϵ_∂ρ∂σ 6
#     @check_size func ρ ∂²ϵ_∂σ² 6
#
#     ccall((:xc_gga_fxc, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}),
#           func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ∂²ϵ_∂ρ²,
#           ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ² )
#     GGASecondDerivative(∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ²)
# end
# function second_energy_derivative{Ρ <: Quantity, ∇Σ <: Quantity}(
#                     func::AbstractLibXCFunctional, ρ::DenseArray{Ρ}, σ::DenseArray{∇Σ})
#     second_energy_derivative!(func,
#                               Units.conversion(Units.ρ, ρ),
#                               Units.conversion(Units.σ, σ),
#                               similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
#                               similar(ρ, Units.∂²ϵ_∂ρ∂σ{Cdouble}, output_size(func, ρ, 6)),
#                               similar(ρ, Units.∂²ϵ_∂σ²{Cdouble}, output_size(func, ρ, 6)))
# end
# function second_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
#                                   σ::DenseArray{Cdouble})
#     second_energy_derivative!(func, ρ, σ,
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
# end
#
# """
#     $(SIGNATURES)
#
# Third derivatives of GGA energy w.r.t. ρ and σ=|σ|. The dimensionality of the arrays is as
# follows:
#
# |GGA        | unpolarized | polarized                 |
# |-----------|-------------|---------------------------|
# |ρ          | any         | `(2, ...)`                |
# |σ          | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
# |∂³ϵ/∂ρ³    | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
# |∂³ϵ/∂ρ²∂σ  | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
# |∂³ϵ/∂ρ∂σ²  | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |
# |∂³ϵ/∂σ³    | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
# """
# function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
#                                   ρ::DenseArray{Units.ρ{Cdouble}},
#                                   σ::DenseArray{Units.σ{Cdouble}},
#                                   ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}},
#                                   ∂³ϵ_∂ρ²∂σ::DenseArray{Units.∂³ϵ_∂ρ²∂σ{Cdouble}},
#                                   ∂³ϵ_∂ρ∂σ²::DenseArray{Units.∂³ϵ_∂ρ∂σ²{Cdouble}},
#                                   ∂³ϵ_∂σ³::DenseArray{Units.∂³ϵ_∂σ³{Cdouble}})
#     third_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂³ϵ_∂ρ³),
#                              reinterpret(Cdouble, ∂³ϵ_∂ρ²∂σ),
#                              reinterpret(Cdouble, ∂³ϵ_∂ρ∂σ²),
#                              reinterpret(Cdouble, ∂³ϵ_∂σ³))
#     GGAThirdDerivative(∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
# end
# function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
#                                   ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
#                                   ∂³ϵ_∂ρ³::DenseArray{Cdouble},
#                                   ∂³ϵ_∂ρ²∂σ::DenseArray{Cdouble},
#                                   ∂³ϵ_∂ρ∂σ²::DenseArray{Cdouble},
#                                   ∂³ϵ_∂σ³::DenseArray{Cdouble})
#     @check_functional func gga
#     @check_availability func kxc
#     @check_size func ρ σ 3
#     @check_size func ρ ∂³ϵ_∂ρ³ 4
#     @check_size func ρ ∂³ϵ_∂ρ²∂σ 9
#     @check_size func ρ ∂³ϵ_∂ρ∂σ² 12
#     @check_size func ρ ∂³ϵ_∂σ³ 10
#
#     ccall((:xc_gga_fxc, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}, Ptr{Cdouble}),
#           func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ∂³ϵ_∂ρ³,
#           ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
#     GGAThirdDerivative(∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ², ∂³ϵ_∂σ³)
# end
# function third_energy_derivative{Ρ <: Quantity, ∇Σ <: Quantity}(
#                 func::AbstractLibXCFunctional, ρ::DenseArray{Ρ}, σ::DenseArray{∇Σ})
#     third_energy_derivative!(func,
#                              Units.conversion(Units.ρ, ρ), Units.conversion(Units.σ, σ),
#                              similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)),
#                              similar(ρ, Units.∂³ϵ_∂ρ²∂σ{Cdouble}, output_size(func, ρ, 9)),
#                              similar(ρ, Units.∂³ϵ_∂ρ∂σ²{Cdouble}, output_size(func, ρ, 12)),
#                              similar(ρ, Units.∂³ϵ_∂σ³{Cdouble}, output_size(func, ρ, 10)))
# end
# function third_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
#                                  σ::DenseArray{Cdouble})
#     second_energy_derivative!(func, ρ, σ,
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 4)),
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 9)),
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 12)),
#                               similar(ρ, eltype(ρ), output_size(func, ρ, 10)))
# end
#
# """
#     $(SIGNATURES)
#
# GGA energy and potential
# """
# function energy_and_potential!(func::AbstractLibXCFunctional{Cdouble},
#                                ρ::DenseArray{Units.ρ{Cdouble}},
#                                σ::DenseArray{Units.σ{Cdouble}},
#                                ϵ::DenseArray{Units.ϵ{Cdouble}},
#                                ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
#                                ∂ϵ_∂σ::DenseArray{Units.∂ϵ_∂σ{Cdouble}})
#     energy_and_potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                           reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
#                           reinterpret(Cdouble, ∂ϵ_∂σ))
#     GGAEnergyAndPotential(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
# end
# function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
#                                σ::DenseArray{Cdouble}, ϵ::DenseArray{Cdouble},
#                                ∂ϵ_∂ρ::DenseArray{Cdouble}, ∂ϵ_∂σ::DenseArray{Cdouble})
#     @check_functional func gga
#     @check_availability func vxc
#     @check_size func ρ σ 3
#     @check_size func ρ ϵ 1
#     @check_size func ρ ∂ϵ_∂ρ 2
#     @check_size func ρ ∂ϵ_∂σ 3
#
#     ccall((:xc_gga_exc_vxc, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}),
#           func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
#     GGAEnergyAndPotential(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ)
# end
# function energy_and_potential(func::AbstractLibXCFunctional, ρ::DenseArray{Units.ρ{Cdouble}},
#                 σ::DenseArray{Units.σ{Cdouble}})
#     energy_and_potential!(func, ρ, σ, similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
#                           similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
#                           similar(ρ, Units.∂ϵ_∂σ{Cdouble}, output_size(func, ρ, 3)))
# end
# function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
#                               ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble})
#     energy_and_potential!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
#                           similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
# end
#
#
# """
#     $(SIGNATURES)
#
# GGA energy, first, second, and third derivatives, in place. Arrays for the first, second,
# and third derivatives are optional. They should be given only if available for that
# particular functional. When requesting higher derivatives, arrays to store the lower
# derivatives should also be given.
# """
# function gga!{T <: DenseArray{Cdouble}}(func::AbstractLibXCFunctional{Cdouble},
#                                         ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
#                                         ϵ::DenseArray{Cdouble}, outputs::Vararg{T})
#     if length(outputs) == 0
#         return AllGGA(energy!(func, ρ, σ, ϵ), [], [], [], [], [], [], [], [])
#     elseif length(outputs) == 2
#         result = energy_and_potential!(func, ρ, σ, ϵ, outputs...)
#         return AllGGA(result[1], result[2], [], [], [], [], [], [], [])
#     end
#
#     @check_functional func gga
#     @check_availability func exc
#     @check_availability func vxc
#     length(outputs) ∉ (5, 9) && throw(ArgumentError("Incorrect number of arguments"))
#     length(outputs) ∈ (5, 9) && @check_availability func fxc
#     length(outputs) == 9 && @check_availability func kxc
#     @check_size func ρ σ 3
#     @check_size func ρ ϵ 1
#     @check_size func ρ outputs[1] 2
#     @check_size func ρ outputs[2] 3
#     @check_size func ρ outputs[3] 3
#     @check_size func ρ outputs[4] 6
#     @check_size func ρ outputs[5] 6
#     if length(outputs) == 9
#         @check_size func ρ outputs[6] 4
#         @check_size func ρ outputs[7] 9
#         @check_size func ρ outputs[8] 12
#         @check_size func ρ outputs[9] 10
#     end
#
#     args = tuple(outputs..., (C_NULL for i in 1:(9 - length(outputs)))...)
#
#     ccall((:xc_gga, libxc), Void,
#           (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
#            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
#           func.c_ptr, length(ρ) / convert(Int64, spin(func)),
#           ρ, σ, ϵ, args[1], args[2], args[3], args[4], args[5], args[6], args[7],
#           args[8], args[9])
#
#     if length(outputs) == 5
#         AllGGA(ϵ, outputs..., [], [], [], [])
#     else
#         AllGGA(ϵ, outputs...)
#     end
# end
# function gga!(func::AbstractLibXCFunctional{Cdouble},
#               ρ::DenseArray{Units.ρ{Cdouble}},
#               σ::DenseArray{Units.σ{Cdouble}},
#               ϵ::DenseArray{Units.ϵ{Cdouble}},
#               ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
#               ∂ϵ_∂σ::DenseArray{Units.∂ϵ_∂σ{Cdouble}},
#               ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
#               ∂²ϵ_∂ρ∂σ::DenseArray{Units.∂²ϵ_∂ρ∂σ{Cdouble}},
#               ∂²ϵ_∂σ²::DenseArray{Units.∂²ϵ_∂σ²{Cdouble}},
#               ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}},
#               ∂³ϵ_∂ρ²∂σ::DenseArray{Units.∂³ϵ_∂ρ²∂σ{Cdouble}},
#               ∂³ϵ_∂ρ∂σ²::DenseArray{Units.∂³ϵ_∂ρ∂σ²{Cdouble}},
#               ∂³ϵ_∂σ³::DenseArray{Units.∂³ϵ_∂σ³{Cdouble}})
#     result = gga!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                   reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
#                   reinterpret(Cdouble, ∂ϵ_∂σ), reinterpret(Cdouble, ∂²ϵ_∂ρ²),
#                   reinterpret(Cdouble, ∂²ϵ_∂ρ∂σ), reinterpret(Cdouble, ∂²ϵ_∂σ²),
#                   reinterpret(Cdouble, ∂³ϵ_∂ρ³), reinterpret(Cdouble, ∂³ϵ_∂ρ²∂σ),
#                   reinterpret(Cdouble, ∂³ϵ_∂ρ∂σ²), reinterpret(Cdouble, ∂³ϵ_∂σ³))
#     AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ², ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂σ, ∂³ϵ_∂ρ∂σ²,
#            ∂³ϵ_∂σ³)
# end
# function gga!(func::AbstractLibXCFunctional{Cdouble},
#               ρ::DenseArray{Units.ρ{Cdouble}}, σ::DenseArray{Units.σ{Cdouble}},
#               ϵ::DenseArray{Units.ϵ{Cdouble}},
#               ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
#               ∂ϵ_∂σ::DenseArray{Units.∂ϵ_∂σ{Cdouble}},
#               ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
#               ∂²ϵ_∂ρ∂σ::DenseArray{Units.∂²ϵ_∂ρ∂σ{Cdouble}},
#               ∂²ϵ_∂σ²::DenseArray{Units.∂²ϵ_∂σ²{Cdouble}})
#     result = gga!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                   reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
#                   reinterpret(Cdouble, ∂ϵ_∂σ), reinterpret(Cdouble, ∂²ϵ_∂ρ²),
#                   reinterpret(Cdouble, ∂²ϵ_∂ρ∂σ), reinterpret(Cdouble, ∂²ϵ_∂σ²))
#     AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂σ, ∂²ϵ_∂σ², Units.∂³ϵ_∂ρ³{Cdouble}[],
#            Units.∂³ϵ_∂ρ²∂σ{Cdouble}[], Units.∂³ϵ_∂ρ∂σ²{Cdouble}[],
#            Units.∂³ϵ_∂σ³{Cdouble}[])
# end
# function gga!(func::AbstractLibXCFunctional{Cdouble},
#               ρ::DenseArray{Units.ρ{Cdouble}}, σ::DenseArray{Units.σ{Cdouble}},
#               ϵ::DenseArray{Units.ϵ{Cdouble}},
#               ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
#               ∂ϵ_∂σ::DenseArray{Units.∂ϵ_∂σ{Cdouble}})
#     result = energy_and_potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                                    reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
#                                    reinterpret(Cdouble, ∂ϵ_∂σ))
#     AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂σ, Units.∂²ϵ_∂ρ²{Cdouble}[], Units.∂²ϵ_∂ρ∂σ{Cdouble}[],
#            Units.∂²ϵ_∂σ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[], Units.∂³ϵ_∂ρ²∂σ{Cdouble}[],
#            Units.∂³ϵ_∂ρ∂σ²{Cdouble}[], Units.∂³ϵ_∂σ³{Cdouble}[])
# end
# function gga!(func::AbstractLibXCFunctional{Cdouble},
#               ρ::DenseArray{Units.ρ{Cdouble}}, σ::DenseArray{Units.σ{Cdouble}},
#               ϵ::DenseArray{Units.ϵ{Cdouble}})
#     result = energy!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, σ),
#                      reinterpret(Cdouble, ϵ))
#     AllGGA(ϵ, Units.∂ϵ_∂ρ{Cdouble}[], Units.∂ϵ_∂σ{Cdouble}[], Units.∂²ϵ_∂ρ²{Cdouble}[],
#            Units.∂²ϵ_∂ρ∂σ{Cdouble}[], Units.∂²ϵ_∂σ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[],
#            Units.∂³ϵ_∂ρ²∂σ{Cdouble}[], Units.∂³ϵ_∂ρ∂σ²{Cdouble}[], Units.∂³ϵ_∂σ³{Cdouble}[])
# end
#
#
# """
#     $(SIGNATURES)
#
# Computes the energy and all available derivatives for the given functional
# """
# function gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}},
#              σ::DenseArray{Units.σ{Cdouble}})
#     @check_functional func gga
#     const f = flags(func)
#     if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
#              similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
#              similar(ρ, Units.∂ϵ_∂σ{Cdouble}, output_size(func, ρ, 3)),
#              similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
#              similar(ρ, Units.∂²ϵ_∂ρ∂σ{Cdouble}, output_size(func, ρ, 6)),
#              similar(ρ, Units.∂²ϵ_∂σ²{Cdouble}, output_size(func, ρ, 6)),
#              similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)),
#              similar(ρ, Units.∂³ϵ_∂ρ²∂σ{Cdouble}, output_size(func, ρ, 9)),
#              similar(ρ, Units.∂³ϵ_∂ρ∂σ²{Cdouble}, output_size(func, ρ, 10)),
#              similar(ρ, Units.∂³ϵ_∂σ³{Cdouble}, output_size(func, ρ, 12)))
#     elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
#              similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
#              similar(ρ, Units.∂ϵ_∂σ{Cdouble}, output_size(func, ρ, 3)),
#              similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
#              similar(ρ, Units.∂²ϵ_∂ρ∂σ{Cdouble}, output_size(func, ρ, 6)),
#              similar(ρ, Units.∂²ϵ_∂σ²{Cdouble}, output_size(func, ρ, 6)))
#     elseif Constants.exc ∈ f && Constants.vxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
#              similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
#              similar(ρ, Units.∂ϵ_∂σ{Cdouble}, output_size(func, ρ, 3)))
#     elseif Constants.exc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
#     else
#         throw(ArgumentError("Not sure what this functional can do"))
#     end
# end
# function gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
#              σ::DenseArray{Cdouble})
#     @check_functional func gga
#
#     const f = flags(func)
#     if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
#              similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 4)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 9)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 12)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 10)))
#     elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
#              similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
#              similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
#     elseif Constants.exc ∈ f && Constants.vxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
#              similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
#     elseif Constants.exc ∈ f && Constants.vxc ∈ f
#         gga!(func, ρ, σ,
#              similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
#              similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
#     elseif Constants.exc ∈ f
#         gga!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
#     else
#         throw(ArgumentError("Not sure what this functional can do"))
#     end
# end
