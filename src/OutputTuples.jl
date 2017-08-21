module OutputTuples
export AllLDA, LDAEnergyAndPotential, GGAPotential, GGAEnergyAndPotential,
       GGASecondDerivative, GGAThirdDerivative, AllGGA, LDATuple, GGATuple
using ..Dispatch
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree
using Base.Iterators: zip
using Unitful

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

""" Energy and potential from LDA """
struct GGAPotential{T0, T1} <: NamedTuple
    ∂ϵ_∂ρ::T0
    ∂ϵ_∂σ::T1
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
    ∂³ϵ_∂ρ∂σ²::T3
    ∂³ϵ_∂σ³::T3
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
Base.start(::NamedTuple) = 1
Base.next(iter::NamedTuple, state::Integer) = iter[state], state + 1
Base.done(iter::NamedTuple, state::Integer) = state > length(iter)

simplify_units(io::IO, a::Any) = show(io, a)
simplify_units{T, D, U}(io::IO, a::AbstractArray{Quantity{T, D, U}}) = begin
    print(io, reinterpret(T, a))
    print(io, "u\"$(U())\"")
end

Base.show(io::IO, t::NamedTuple) = begin
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
        $name(args::Vararg{DD.AxisArrays.All}) = begin
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
