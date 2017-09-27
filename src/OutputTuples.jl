module OutputTuples
export AllLDA, LDAEnergyAndPotential, GGAPotential, GGAEnergyAndPotential,
       GGASecondDerivative, GGAThirdDerivative, AllGGA, LDATuple, GGATuple,
       SLDAEnergyAndPotential, SLDASecondDerivative, SLDAThirdDerivative,
       GGASecondDerivatives, GGAPotentials, GGAThirdDerivatives
using ..Dispatch
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree
using Base.Iterators: zip
using Unitful
using AxisArrays
using StaticArrays

macro lintpragma(s) end

abstract type NamedTuple end

""" All outputs from LDA """
struct AllLDA{T0, T1, T2, T3} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
    ∂²ϵ_∂ρ²::T2
    ∂³ϵ_∂ρ³::T3
end

""" Energy and potential from LDA """
struct LDAEnergyAndPotential{T0, T1} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
end

""" Energy and potential from spin-LDA """
struct SLDAEnergyAndPotential{T0, T1} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρα::T1
    ∂ϵ_∂ρβ::T1
end

""" Second energy derivative for spin-LDA """
struct SLDASecondDerivative{T} <: FieldVector{3, T}
    ∂²ϵ_∂ρα²::T
    ∂²ϵ_∂ρα∂ρβ::T
    ∂²ϵ_∂ρβ²::T
end

""" Third energy derivative for spin-LDA """
struct SLDAThirdDerivative{T} <: FieldVector{4, T}
    ∂³ϵ_∂ρα³::T
    ∂³ϵ_∂ρα²∂ρβ::T
    ∂³ϵ_∂ρα∂ρβ²::T
    ∂³ϵ_∂ρβ³::T
end

""" Energy and potential from LDA """
struct GGAPotential{T0, T1} <: NamedTuple
    ∂ϵ_∂ρ::T0
    ∂ϵ_∂σ::T1
end

struct GGAPotentials{T0, T1} <: NamedTuple
    ∂ϵ_∂ρα::T0
    ∂ϵ_∂ρβ::T0
    ∂ϵ_∂σαα::T1
    ∂ϵ_∂σαβ::T1
    ∂ϵ_∂σββ::T1
end

struct GGASecondDerivatives{T0, T1, T2} <: NamedTuple
    ∂²ϵ_∂ρα²::T0
    ∂²ϵ_∂ρα∂ρβ::T0
    ∂²ϵ_∂ρβ²::T0
    ∂²ϵ_∂ρα∂σαα::T1
    ∂²ϵ_∂ρα∂σαβ::T1
    ∂²ϵ_∂ρα∂σββ::T1
    ∂²ϵ_∂ρβ∂σαα::T1
    ∂²ϵ_∂ρβ∂σαβ::T1
    ∂²ϵ_∂ρβ∂σββ::T1
    ∂²ϵ_∂σαα²::T2
    ∂²ϵ_∂σαα∂σαβ::T2
    ∂²ϵ_∂σαα∂σββ::T2
    ∂²ϵ_∂σαβ²::T2
    ∂²ϵ_∂σαβ∂σββ::T2
    ∂²ϵ_∂σββ²::T2
end

""" Second derivative from GGA

Include the second derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
struct GGASecondDerivative{T0, T1, T2} <: NamedTuple
    ∂²ϵ_∂ρ²::T0
    ∂²ϵ_∂ρ∂σ::T1
    ∂²ϵ_∂σ²::T2
end

""" Third derivative from GGA

Include the third derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
struct GGAThirdDerivative{T0, T1, T2, T3} <: NamedTuple
    ∂³ϵ_∂ρ³::T0
    ∂³ϵ_∂ρ²∂σ::T1
    ∂³ϵ_∂ρ∂σ²::T2
    ∂³ϵ_∂σ³::T3
end
struct GGAThirdDerivatives{T0, T1, T2, T3} <: NamedTuple
    ∂³ϵ_∂ρα³::T0
    ∂³ϵ_∂ρα²∂ρβ::T0
    ∂³ϵ_∂ρα∂ρβ²::T0
    ∂³ϵ_∂ρβ³::T0

    ∂³ϵ_∂ρα²∂σαα::T1
    ∂³ϵ_∂ρα²∂σαβ::T1
    ∂³ϵ_∂ρα²∂σββ::T1
    ∂³ϵ_∂ρα∂ρβ∂σαα::T1
    ∂³ϵ_∂ρα∂ρβ∂σαβ::T1
    ∂³ϵ_∂ρα∂ρβ∂σββ::T1
    ∂³ϵ_∂ρβ²∂σαα::T1
    ∂³ϵ_∂ρβ²∂σαβ::T1
    ∂³ϵ_∂ρβ²∂σββ::T1

    ∂³ϵ_∂ρα∂σαα²::T2
    ∂³ϵ_∂ρα∂σαα∂σαβ::T2
    ∂³ϵ_∂ρα∂σαα∂σββ::T2
    ∂³ϵ_∂ρα∂σαβ²::T2
    ∂³ϵ_∂ρα∂σαβ∂σββ::T2
    ∂³ϵ_∂ρα∂σββ²::T2
    ∂³ϵ_∂ρβ∂σαα²::T2
    ∂³ϵ_∂ρβ∂σαα∂σαβ::T2
    ∂³ϵ_∂ρβ∂σαα∂σββ::T2
    ∂³ϵ_∂ρβ∂σαβ²::T2
    ∂³ϵ_∂ρβ∂σαβ∂σββ::T2
    ∂³ϵ_∂ρβ∂σββ²::T2

    ∂³ϵ_∂σαα³::T3
    ∂³ϵ_∂σαα²∂σαβ::T3
    ∂³ϵ_∂σαα²∂σββ::T3
    ∂³ϵ_∂σαα∂σαβ²::T3
    ∂³ϵ_∂σαα∂σαβ∂σββ::T3
    ∂³ϵ_∂σαα∂σββ²::T3
    ∂³ϵ_∂σαβ³::T3
    ∂³ϵ_∂σαβ²∂σββ::T3
    ∂³ϵ_∂σαβ∂σββ²::T3
    ∂³ϵ_∂σββ³::T3
end

""" Holds GGA energy and first derivatives """
struct GGAEnergyAndPotential{T0, T1, T2} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
    ∂ϵ_∂σ::T2
end

""" All outputs from LDA """
struct AllGGA{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
    ∂ϵ_∂σ::T2
    ∂²ϵ_∂ρ²::T3
    ∂²ϵ_∂ρ∂σ::T4
    ∂²ϵ_∂σ²::T5
    ∂³ϵ_∂ρ³::T6
    ∂³ϵ_∂ρ²∂σ::T7
    ∂³ϵ_∂ρ∂σ²::T8
    ∂³ϵ_∂σ³::T9
end

Base.length(a::NamedTuple) = length(fieldnames(typeof(a)))
Base.getindex(a::NamedTuple, index::Integer) =
    getfield(a, fieldname(typeof(a), index))
Base.getindex(a::NamedTuple, sequence::Union{OrdinalRange, Vector{<:Integer}}) =
    tuple((getfield(a, fieldname(typeof(a), i)) for i in sequence)...)
Base.start(::NamedTuple) = 1
Base.next(iter::NamedTuple, state::Integer) = iter[state], state + 1
Base.done(iter::NamedTuple, state::Integer) = state > length(iter)
Base.endof(a::NamedTuple) = length(a)

simplify_units(io::IO, a::Any) = show(io, a)
simplify_units(io::IO, a::AxisArray) = simplify_units(io, a.data)
simplify_units(io::IO, a::AbstractArray{<: Quantity}) = begin
    print(io, reinterpret(eltype(one(eltype(a))), a))
    print(io, "u\"$(unit(eltype(a)))\"")
end

Base.show(io::IO, t::Union{NamedTuple, SLDASecondDerivative, SLDAThirdDerivative}) = begin
    print(io, "(")
    isfirst = true
    for (k,v) in zip(fieldnames(t), t)
        !isfirst && print(io, ", ")
        print(io, k, " = "); simplify_units(io, v)
        isfirst = false
    end
    print(io, ")")
end

for (cons, name) in (:AllLDA => :LDATuple, :AllGGA => :GGATuple)
    @eval begin
        $name(args::Vararg{DD.Arrays.All}) = begin
            @lintpragma("Ignore unused i")
            concretize = i -> i{Int64}
            fielddims = dimension.(concretize.(getfield.(DH.Scalars, fieldnames($cons))))

            reordered = Any[nothing for i in 1:length(fieldnames($cons))]
            any(reordered .== 0) && throw("Argument with unexpected Unitful dimension")
            for arg in args
                j = findfirst(fielddims, dimension(eltype(arg)))
                j == 0 && throw("Argument with unexpected Unitful dimension")
                reordered[j] = arg
            end
            $cons{typeof.(reordered)...}(reordered...)
        end
    end
end

end
