module FunctionalMacros
import LibXC: energy, energy!, potential, potential!, energy_and_potential,
              energy_and_potential!, second_energy_derivative, second_energy_derivative!,
              third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Constants
using ..Internals
using ..Internals: AbstractLibXCFunctional, CFuncType, libxc
using ..ArrayManips: from_libxc_array!, to_libxc_array, valid_array
using ..OutputTuples 

using AxisArrays
using ArgCheck
using Unitful

using DFTShims: ColinearSpinFirst, SpinDegenerate, Dispatch, is_spin_polarized, components, 
                SpinCategory, SpinAware
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

_requested_outputs(o::Vararg{DD.Arrays.All}) = begin
    @lintpragma("Ignore unused o")
    check = (x, y) -> eltype(x) <: y
    result = Constants.FLAGS[]
    eltypes = eltype.(o)
    findfirst(x -> x <: DD.Scalars.ϵ, eltypes) ≠ 0 && push!(result, Constants.exc)
    findfirst(x -> x <: DD.Scalars.∂ϵ_∂ρ, eltypes) ≠ 0 && push!(result, Constants.vxc)
    findfirst(x -> x <: DD.Scalars.∂²ϵ_∂ρ², eltypes) ≠ 0 && push!(result, Constants.fxc)
    findfirst(x -> x <: DD.Scalars.∂³ϵ_∂ρ³, eltypes) ≠ 0 && push!(result, Constants.kxc)
    unique(result)
end

""" Constructs mutating functions taking AxisArray arguments """
_mutating_wrapper_functionals(name::Symbol, dfttype::Symbol, output_type::Symbol,
                              outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
   
    allargs = dfttype == :lda ? (:ρ, outputs...): (:ρ, :σ, outputs...)

    # function arguments
    axis_args = (arg -> :($arg::DD.AxisArrays.$arg)).(allargs)
    array_args = (arg -> :($arg::DD.DenseArrays.$arg{Cdouble})).(allargs)

    # converts backt to original axisarrays
    name! = Symbol(name, :!)
    # input sanity checks
    arg_checks = x -> :(valid_array(ρ, $x))
    arg_checks_array = x -> :(valid_array(ρ, $x, Spin))
    # converts to type understood by C LibXC
    to_libxc_types = x -> :($(Symbol(:c, x)) = to_libxc_array(ρ, $x))
    to_libxc_array_types = x -> :($(Symbol(:c, x)) = to_libxc_array($x))
    to_dense_types = x -> :(reinterpret(Cdouble, $(Symbol(:c, x)).data))
    to_c_types = x -> :(reinterpret(Cdouble, $(Symbol(:c, x))))
    from_libxc_types = x -> :(from_libxc_array!($x, $(Symbol(:c, x))))

    quote
        $(esc(name!))(func::AbstractLibXCFunctional, $(axis_args...)) = begin
            $(arg_checks.(allargs)...)
            if family(func) ≠ getfield(Constants, $(QuoteNode(dfttype)))
                throw(ArgumentError("Incorrect functional type $($(QuoteNode(dfttype)))"))
            end
            reqout = _requested_outputs($(outputs...))
            if !(reqout ⊆ flags(func))
                throw(ArgumentError("This functional does not implement all of $reqout"))
            end

            $(to_libxc_types.(allargs)...)
            $(esc(name!))(func, $(to_dense_types.(allargs)...))
            $(from_libxc_types.(outputs)...)
            $(esc(output_type))($(outputs...))
        end

        $(esc(name!))(func::AbstractLibXCFunctional, $(array_args...)) = begin
            const Spin = SpinCategory(func)
            $(arg_checks_array.(allargs)...)
            if family(func) ≠ getfield(Constants, $(QuoteNode(dfttype)))
                throw(ArgumentError("Incorrect functional type $($(QuoteNode(dfttype)))"))
            end
            reqout = _requested_outputs($(outputs...))
            if !(reqout ⊆ flags(func))
                throw(ArgumentError("This functional does not implement all of $reqout"))
            end

            $(to_libxc_array_types.(allargs)...)
            $(esc(name!))(func, $(to_c_types.(allargs)...))
            $(esc(output_type))($(from_libxc_types.(outputs)...))
        end
    end
end

""" Constructs non-mutating functions taking AxisArray arguments """
_nonmutating_wrapper_functionals(name::Symbol, dfttype::Symbol,
                                 outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
    @lintpragma("Ignore unused outputs")

    name! = Symbol(name, :!)
    inputs = dfttype == :lda ? (:ρ, ) : (:ρ, :σ)
    args_axis = x -> :($x::DD.AxisArrays.$x)
    args_array = x -> :($x::DD.Arrays.$x)
    similars = x -> begin
        if x == :ϵ
            :(similar(ρ, DH.Scalars.$x{Cdouble}, SpinDegenerate()))
        else
            :(similar(ρ, DH.Scalars.$x{Cdouble}, SpinAware()))
        end
    end
    similars_nospin = x -> begin
        @lintpragma("Ignore unused x")
        :(similar(ρ, DH.Scalars.$x{Cdouble}))
    end
    similars_spin = x -> begin
        if x == :ϵ
            :(similar(ρ, DH.Scalars.$x{Cdouble},
                      ifelse(ndims(ρ) == 1, (1, ), Base.tail(size(ρ)))))
        else
            :(similar(ρ, DH.Scalars.$x{Cdouble},
                      length(components(DH.Scalars.$x, ColinearSpinFirst())),
                      Base.tail(size(ρ))...))
        end
    end

    quote
        $(esc(name))(func::AbstractLibXCFunctional, $(args_axis.(inputs)...)) =
            $(esc(name!))(func, $(inputs...), $(similars.(outputs)...))
        $(esc(name))(name::Symbol, $(args_axis.(inputs)...)) =
            $(esc(name))(XCFunctional(name, is_spin_polarized(ρ)), $(inputs...))
        $(esc(name))(func::AbstractLibXCFunctional, $(args_array.(inputs)...)) = begin
            if spin(func) == Constants.polarized
                $(esc(name!))(func, $(inputs...), $(similars_spin.(outputs)...))
            else
                $(esc(name!))(func, $(inputs...), $(similars_nospin.(outputs)...))
            end
        end
    end
end

""" Constructs functions which actually call the C implementations """
_c_wrapper_functionals(name::Symbol, dfttype::Symbol, cname::Symbol,
                       outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
    allargs = dfttype == :lda ? (:ρ, outputs...): (:ρ, :σ, outputs...)
    name! = Symbol(name, :!)
    nonnull = collect(Iterators.filter(x -> x ≠ :C_NULL, allargs))
    dargs = (arg -> :($arg::DenseArray{Cdouble})).(nonnull)
    quote
        $(esc(name!))(func::AbstractLibXCFunctional, $(dargs...)) = begin
            @lintpragma("Ignore use of undeclared variable ρ")
            @lintpragma("Ignore unused u")
            @lintpragma("Ignore use of undeclared variable ccall")
            ccall(($(QuoteNode(cname)), $libxc), Void,
                  (Ptr{CFuncType}, Cint, $([:(Ptr{Cdouble}) for u in allargs]...)),
                  func.c_ptr, length(ρ) / (spin(func) == Constants.polarized ? 2: 1),
                  $(allargs...))
            tuple($(nonnull...))
        end
    end
end

""" Constructs scalar functions """
macro _scalar_functional(name::Symbol, _inputs::Expr, _outputs::Expr, output_type::Symbol)
    unit_args = arg -> :($(arg.args[1])::DD.Scalars.$(arg.args[2]))
    args = arg -> arg.args[1]
    DD2Float = arg ->
        :(ustrip(convert(DH.Scalars.$(arg.args[2]){Cdouble}, $(arg.args[1]))))
    Float2DD = arg -> :(DH.Scalars.$(arg[2]){T}(result[$(arg[1])]))
    types = arg -> :(eltype(one($(arg.args[1]))))
    inputs = tuple(_inputs.args...)
    outputs = tuple(_outputs.args...)
    polarization = countnz(x.args[2] == :ρ for x in inputs) == 1 ?
                    Constants.unpolarized:
                    Constants.polarized
    dfttype = countnz(x.args[2] == :σ for x in inputs) ≥ 1 ?
                    Constants.gga:
                    Constants.lda
    private_name = Symbol(:unsafe_, name)
    checks = name ∈ (:xc_lda_exc, :xc_gga_exc) ? (Constants.exc, ):
             name ∈ (:xc_lda_vxc, :xc_gga_vxc) ? (Constants.vxc, ):
             name ∈ (:xc_lda_fxc, :xc_gga_fxc) ? (Constants.fxc, ):
             name ∈ (:xc_lda_kxc, :xc_gga_kxc) ? (Constants.kxc, ):
                 (Constants.exc, Constants.vxc)
    # because otherwise, inferrence fails
    first_few = Base.first(Base.IteratorsMD.split(inputs, Val{2}))
    quote
        $(esc(name))(name::Symbol, $(unit_args.(inputs)...)) =
        $(esc(name))(XCFunctional(name, $polarization), $(args.(inputs)...))
        $(esc(name))(func::AbstractLibXCFunctional, $(unit_args.(inputs)...)) = begin
            @argcheck family(func) == $dfttype
            @argcheck spin(func) == $polarization
            @argcheck $checks ⊆ flags(func)
            result = $(esc(private_name))(func, $(DD2Float.(inputs)...))
            T = promote_type(Cdouble, $(types.(first_few)...))
            $output_type($(Float2DD.(collect(enumerate(outputs)))...))
        end
    end
end

macro _all_wrapper_functionals(func, ftype, cname, output_type, outputs...)
    quote
        $(_mutating_wrapper_functionals(func, ftype, output_type, outputs))
        $(_nonmutating_wrapper_functionals(func, ftype, outputs))
        $(_c_wrapper_functionals(func, ftype, cname, outputs))
    end
end

macro _wrapper_functionals(n::Int, func, cname, output_type, outputs...)
    @lintpragma("Ignore unused i")
    o = (outputs..., (:C_NULL for i in (length(outputs) + 1):n)...)
    quote
        $(_mutating_wrapper_functionals(func, func, output_type, outputs))
        $(_c_wrapper_functionals(func, func, cname, o))
    end
end

end
