"""
LibXC compatible dimensions for dispatch

Using these types might result in conversion to appropriate floating point and units, e.g.
copies and such.
"""
module DUnits
using Unitful
using UnitfulHartree

macro _dim_helper(name, quant)
    q = eval(quant)
    dims = typeof(dimension(q))
    units = typeof(unit(q))
    esc(quote
            Unitful.Compat.@compat ($name){T,U} = Unitful.Quantity{T, $dims, U}
            $name(val::Number, units::Unitful.Units=$(units())) =
                $name{typeof(val), typeof(units)}(val)
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
