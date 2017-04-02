""" LibXC compatible units """
module Units
using DocStringExtensions
using Unitful
using UnitfulHartree

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


macro _dim_helper(name, quant)
    q = eval(quant)
    dims = typeof(dimension(q))
    units = typeof(unit(q))
    esc(quote
        Unitful.Compat.@compat ($name){T} = Unitful.Quantity{T, $dims, $units}
        $name(val::Number) = $name{typeof(val)}(val)
        Unitful.unit(::Type{$name}) = $units()
        end)
end

@_dim_helper ρ          1u"𝐞*a₀^-3"
@_dim_helper ∇ρ         1u"𝐞*a₀^-4"
@_dim_helper ϵ          1u"Eₕ/𝐞"
@_dim_helper ∂ϵ_∂ρ      1u"Eₕ*𝐞^-2*a₀^3"
@_dim_helper ∂ϵ_∂∇ρ     1u"Eₕ*𝐞^-2*a₀^4"
@_dim_helper ∂²ϵ_∂ρ²    1u"Eₕ*𝐞^-3*a₀^6"
@_dim_helper ∂²ϵ_∂∇ρ²   1u"Eₕ*𝐞^-3*a₀^8"
@_dim_helper ∂²ϵ_∂ρ∂∇ρ  1u"Eₕ*𝐞^-3*a₀^7"
@_dim_helper ∂³ϵ_∂ρ³    1u"Eₕ*𝐞^-4*a₀^9"
@_dim_helper ∂³ϵ_∂∇ρ³   1u"Eₕ*𝐞^-4*a₀^12"
@_dim_helper ∂³ϵ_∂ρ²∂∇ρ 1u"Eₕ*𝐞^-4*a₀^10"
@_dim_helper ∂³ϵ_∂ρ∂∇ρ² 1u"Eₕ*𝐞^-4*a₀^11"

end
