""" LibXC compatible units

Helps with dispatching over different units to hartree units.
"""
module Units
using DocStringExtensions
using Unitful
using UnitfulHartree
using LibXC.DFTUnits: 𝐞
using LibXC: DFTUnits

"""
    $(SIGNATURES)

Convert to given units and Cdouble. This function is meant to simplify the process of
creating a valid array for LibXC.
"""
function conversion{Q <: Unitful.Quantity}(u::Unitful.Units, input::Array{Q})
    factor = uconvert(u, Q(1))
    result = similar(input,
                     Quantity{Cdouble, typeof(dimension(factor)), typeof(unit(factor))})
    copy!(ustrip(result), ustrip(input) * ustrip(factor))
    result
end

"""
    $(SIGNATURES)

Assumes given units and converts to Cdouble. This function is meant to simplify the process
of creating a valid array for LibXC.
"""
function conversion{N <: Number}(u::Unitful.Units, input::Array{N})
    result = similar(input, typeof(one(Cdouble)u))
    copy!(ustrip(result), input)
    result
end

conversion{Q <: Unitful.Quantity}(::Type{Q}, input::Array) = conversion(unit(Q(1)), input)

macro _dim_helper(name, units)
    dims = typeof(dimension(eval(units)))
    tunits = typeof(eval(units))
    esc(quote
        Unitful.Compat.@compat ($name){T} = Unitful.Quantity{T, $dims, $tunits}
        $name(val::Number) = $name{typeof(val)}(val)
        Unitful.unit(::Type{$name}) = $tunits()
        end)
end

@_dim_helper ρ          DFTUnits.ρ
@_dim_helper ∇ρ         DFTUnits.∇ρ
@_dim_helper ϵ          DFTUnits.ϵ
@_dim_helper ∂ϵ_∂ρ      DFTUnits.∂ϵ_∂ρ
@_dim_helper ∂ϵ_∂∇ρ     DFTUnits.∂ϵ_∂∇ρ
@_dim_helper ∂²ϵ_∂ρ²    DFTUnits.∂²ϵ_∂ρ²
@_dim_helper ∂²ϵ_∂∇ρ²   DFTUnits.∂²ϵ_∂∇ρ²
@_dim_helper ∂²ϵ_∂ρ∂∇ρ  DFTUnits.∂²ϵ_∂ρ∂∇ρ
@_dim_helper ∂³ϵ_∂ρ³    DFTUnits.∂³ϵ_∂ρ³
@_dim_helper ∂³ϵ_∂∇ρ³   DFTUnits.∂³ϵ_∂∇ρ³
@_dim_helper ∂³ϵ_∂ρ²∂∇ρ DFTUnits.∂³ϵ_∂ρ²∂∇ρ
@_dim_helper ∂³ϵ_∂ρ∂∇ρ² DFTUnits.∂³ϵ_∂ρ∂∇ρ²

end
