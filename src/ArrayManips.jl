module ArrayManips

using AxisArrays, ArgCheck, Unitful


using DFTShims: ColinearSpinFirst, SpinDegenerate, Dispatch, is_spin_polarized, components, 
                SpinCategory, SpinAware, ColinearSpin, standard_name, ColinearSpin
using DFTShims.Traits: hartree_concretize_type
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

"""
Checks ρ and other array are compatible

Note that the energy ϵ is not spin-polarized, whereas the potential is. Both share the same
physical units. This means we have to make a special case of ϵ.
"""
valid_array(ρ::DD.AxisArrays.ρ, other::DD.AxisArrays.ϵ) = begin
    @lintpragma("Ignore use of undeclared variable u")
    varname = standard_name(other)
    if is_spin_polarized(ρ)
        @argcheck ndims(ρ) == (ndims(other) + 1) "Dimensions of ρ and $varname do not match"
    else
        @argcheck ndims(ρ) == ndims(other) "Dimensions of ρ and $varname do not match"
    end
    for name in (u for u in axisnames(ρ) if u ≠ :spin)
        @argcheck name ∈ axisnames(other)
        @argcheck size(ρ, Axis{name}) == size(other, Axis{name})
    end
end
valid_array(ρ::DD.AxisArrays.ρ, other::DD.AxisArrays.All) = begin
    @lintpragma("Ignore use of undeclared variable u")
    name = standard_name(other)
    @argcheck(is_spin_polarized(other) == is_spin_polarized(ρ),
              "Spin polarization of ρ and $name do not match")
    @argcheck ndims(ρ) == ndims(other) "Dimensions of ρ and $name do not match"
    if is_spin_polarized(ρ)
        ncomps = length(components(eltype(other), ColinearSpin()))
        @argcheck(size(other, Axis{:spin}) == ncomps,
                  "Length of spin axis of $name is incorrect")
    end
    for name in (u for u in axisnames(ρ) if u ≠ :spin)
        @argcheck name ∈ axisnames(other)
        @argcheck size(ρ, Axis{name}) == size(other, Axis{name})
    end
end
valid_array(ρ::DD.DenseArrays.ρ, other::DD.DenseArrays.All, ::SpinDegenerate) = 
    @argcheck(size(ρ) == size(other),
              "Dimensions of ρ and $(standard_name(other)) do not match")
valid_array(ρ::DD.DenseArrays.ρ, other::DD.DenseArrays.ϵ, ::ColinearSpinFirst) =  begin
    if ndims(ρ) == 1
        @argcheck size(other) == (1, ) "Dimensions of ρ and ϵ do not match"
    else
        @argcheck Base.tail(size(ρ)) == size(other) "Dimensions of ρ and ϵ do not match"
    end
end
valid_array(ρ::DD.DenseArrays.ρ, other::DD.DenseArrays.All, C::ColinearSpinFirst) = begin
    @argcheck(Base.tail(size(ρ)) == Base.tail(size(other)),
              "Dimensions of ρ and $(standard_name(other)) do not match")
    @argcheck(size(other, 1) == length(components(eltype(other), C)),
              "Spin dimensions of $(standard_name(other)) is incorrect")
end

""" True if axes are the same """
_libxc_axes(a::DD.AxisArrays.All, b::DD.AxisArrays.All) = begin
    if is_spin_polarized(a) 
        axes(view(a, Axis{:spin}(1))) == axes(view(b, Axis{:spin}(1))) &&
            SpinCategory(b) isa ColinearSpinFirst
    else
        axes(a) == axes(b)
    end
end
_libxc_axes(a::DD.AxisArrays.ϵ, b::DD.AxisArrays.ϵ) = axes(a) == axes(b)
_libxc_axes(a::DD.AxisArrays.ϵ, b::DD.AxisArrays.All) = _libxc_axes(b, a)
_libxc_axes(a::DD.AxisArrays.All, b::DD.AxisArrays.ϵ) =
    is_spin_polarized(a) ? axes(view(a, Axis{:spin}(1))) == axes(b):
                           axes(a) == axes(b)

""" The spin trait LibXC requires """
_libxc_spin(C::SpinDegenerate) = C
_libxc_spin(::ColinearSpin) = ColinearSpinFirst()

""" Converts input array to C libxc array """
to_libxc_array(ρ::DD.AxisArrays.ρ, array::DD.AxisArrays.All) =  begin
    const T = hartree_concretize_type(eltype(array))
    const Q = Quantity{Cdouble, typeof(dimension(T)), typeof(unit(T))}
    eltype(array) === Q && _libxc_axes(ρ, array) && return array
    copy!(similar(ρ, Q, _libxc_spin(SpinCategory(array))), array)
end

to_libxc_array(array::DH.DenseArrays.All{Cdouble}) = array
to_libxc_array(array::DD.DenseArrays.All) = begin
    const T = hartree_concretize_type(eltype(array))
    const Q = Quantity{Cdouble, typeof(dimension(T)), typeof(unit(T))}
    copy!(similar(array, Q), array)
end

""" Converts back from C libxc array """
from_libxc_array!(output::DD.Arrays.All, carray::AbstractArray) = begin
    carray !== output && copy!(output, carray)
    output
end
end
