module LDA

import LibXC: energy, energy!, potential, potential!, energy_and_potential,
              energy_and_potential!, second_energy_derivative, second_energy_derivative!,
              third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Constants
using ..Internals
using ..Internals: CFuncType, libxc, XCFunctional
using ..OutputTuples
using ..Checks

using AxisArrays
using DocStringExtensions

using DFTShims: SpinCategory, SpinDegenerate, ColinearSpinFirst, Dispatch,
                is_spin_polarized, components
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable ccall")


const TYPES = (:ϵ, :∂ϵ_∂ρ, :∂²ϵ_∂ρ², :∂³ϵ_∂ρ³)
const FUNCTIONS = Dict(:energy => TYPES[1:1], :energy_and_potential => TYPES[1:2],
                       :potential => TYPES[2:2], :second_energy_derivative => TYPES[3:3],
                       :third_energy_derivative => TYPES[4:4], :lda => TYPES)
const OUTPUTS = Dict(:energy => :AxisArray, :energy_and_potential => :LDAEnergyAndPotential,
                     :potential => :AxisArray, :second_energy_derivative => :AxisArray,
                     :third_energy_derivative => :AxisArray, :lda => :LDATuple)
const FUNCNAMES = Dict(:energy => :xc_lda_exc, :energy_and_potential => :xc_lda_exc_vxc,
                       :potential => :xc_lda_vxc, :second_energy_derivative => :xc_lda_fxc,
                       :third_energy_derivative => :xc_lda_kxc, :lda => :xc_lda)
const FUNCTYPES = Dict(:energy => :exc, :energy_and_potential => (:exc, :vxc),
                       :potential => :vxc, :second_energy_derivative => :fxc,
                       :third_energy_derivative => :kxc, :lda => (:exc, :vxc, :fxc, :kxc))


# Argument of functions constructed later
_ddargument(arg::Symbol) = :($arg::DD.AxisArrays.$arg)
_dhargument(arg::Symbol) = :($arg::DH.AxisArrays.$arg{Float64})
_deargument(arg::Symbol) = :($arg::DenseArray{Float64})
# Conversion between AxisArray of different unit/types to something LibXC understands
_conversion(etype::Type, array::AxisArray) = begin
    eltype(array) === etype && return array
    result = similar(array, etype)
    result .= array
    result
end
_conversion(arg::Symbol) = :(_conversion(DH.Scalars.$arg{Float64}, $arg))
# Conversion back from LibXC specific types
_convert_back(arg::Symbol) =
    :(AxisArray(convert(typeof($arg.data), result.$arg.data), axes($arg)))
_convert_back(arg::Tuple{Symbol}) =
    (:(AxisArray(convert(typeof($(arg[1]).data), result.data), axes($(arg[1])))), )
_convert_back(args::Tuple{Symbol, Vararg{Symbol}}) = begin
    @lintpragma("Ignore unused args")
    _convert_back.(args)
end
# Removes AxisArray wrapper before call to no-frill function
_convert_to_array(arg::Symbol) =
    :(reinterpret(typeof(one(eltype($arg))), convert(Array, $arg)))
# Helps create output arrays in unbangued functions (e.g. when calling lda! from lda)
_similar(arg::Symbol) = begin
    @lintpragma("Ignore unused arg")
    :(similar(DH.Scalars.$arg{Float64}, ρ))
end
# Inserts assertions to check correctness of input data during calls
_check_spin(b::Symbol, a::Symbol) = begin
    msg = "Spin of $a and $b do not match"
    :(SpinCategory($a) == SpinCategory($b) || throw(ArgumentError($msg)))
end
_check_axes_nospin(b::Symbol) = begin
    msg = "Axes of ρ and $b do not match"
    :(all(axes(ρ) .== axes($b)) || throw(ArgumentError($msg)))
end
_check_axes_spin(b::Symbol) = begin
    msg = "Axes of ρ and $b do not match"
    spinmsg = "First axis of $b is not spin"
    spinerror = "Axis values of spin axis of $b are incorrect"
    comps = gensym("comps")
    quote
        $comps = components(eltype($b), ColinearSpinFirst())
        if length($comps) > 1
            axisnames($b)[1] ≠ :spin && throw(ArgumentError($spinmsg))
            axisvalues($b)[1] ≠ $comps && throw(ArgumentError($spinerror))
            all(axes(ρ)[2:end] .== axes($b)[2:end]) || throw(ArgumentError($msg))
        end
    end
end

_check_functional(functype::Symbol) = begin
    msg = "Functional is LDA, but input array do not correspond"
    :(family($functype) == Constants.lda || throw(ArgumentError($msg)))
end
_check_availability(out::Symbol, name::Symbol) = begin
    msg = "This functional does not implement $out"
    :(Constants.$out ∈ flags($name) || error($msg))
end

for (func, outputs) in FUNCTIONS
    allargs = (:ρ, outputs...)
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
            $(_check_availability.(FUNCTYPES[func], :func))

            $func!(func, $(_convert_to_array.(allargs)...))

            $otype($(outputs...))
        end

        $func!(func::AbstractLibXCFunctional, $(_deargument.(allargs)...)) = begin
            @lintpragma("Ignore use of undeclared variable ρ")
            @lintpragma("Ignore unused u")
            ccall(($(QuoteNode(FUNCNAMES[func])), $libxc), Void,
                  (Ptr{CFuncType}, Cint, $([:(Ptr{Float64}) for u in allargs]...)),
                  func.c_ptr, length(ρ) / (spin(func) == Constants.polarized ? 2: 1),
                  ρ, $(outputs...))
            $outputs
        end

        $func(func::AbstractLibXCFunctional, ρ::DD.AxisArrays.ρ) =
            $func!(func, ρ, $(_similar.(outputs)...))

        $func(name::Symbol, ρ::DD.AxisArrays.ρ) =
            $func(XCFunctional(name, is_spin_polarized(ρ)), ρ)
    end

    func == :lda && continue

    @eval lda!(func::AbstractLibXCFunctional, $(_ddargument.(allargs)...)) = begin
        result = $func(func, $(allargs...))
        @lintpragma("Ignore use of undeclared variable LDATuple")
        LDATuple(result...)
    end
end

"""
    $(SIGNATURES)

Computes the energy and all available derivatives for the given functional
"""
lda(func::AbstractLibXCFunctional{Float64}, ρ::DH.AxisArrays.ρ{Float64}) = begin
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
