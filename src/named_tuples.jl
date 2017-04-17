abstract NamedTuple

""" All outputs from LDA """
immutable AllLDA{T0, T1, T2, T3} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
    ∂²ϵ_∂ρ²::T2
    ∂³ϵ_∂ρ³::T3
end
Base.length(a::AllLDA) = 4
function Base.getindex(a::AllLDA, index::Integer)
    if index == 1
        a.ϵ
    elseif index == 2
        a.∂ϵ_∂ρ
    elseif index == 3
        a.∂²ϵ_∂ρ²
    elseif index == 4
        a.∂³ϵ_∂ρ³
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" Energy and potential from LDA """
immutable LDAEnergyAndPotential{T0, T1} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
end
Base.length(a::LDAEnergyAndPotential) = 2
function Base.getindex(a::LDAEnergyAndPotential, index::Integer)
    if index == 1
        a.ϵ
    elseif index == 2
        a.∂ϵ_∂ρ
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" Energy and potential from LDA """
immutable GGAPotential{T0, T1} <: NamedTuple
    ∂ϵ_∂ρ::T0
    ∂ϵ_∂σ::T1
end
Base.length(a::GGAPotential) = 2
function Base.getindex(a::GGAPotential, index::Integer)
    if index == 1
        a.∂ϵ_∂ρ
    elseif index == 2
        a.∂ϵ_∂σ
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" Second derivative from GGA

Include the second derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
immutable GGASecondDerivative{T0, T1, T2} <: NamedTuple
    ∂²ϵ_∂ρ²::T0
    ∂²ϵ_∂ρ∂σ::T1
    ∂²ϵ_∂σ²::T2
end
Base.length(a::GGASecondDerivative) = 3
function Base.getindex(a::GGASecondDerivative, index::Integer)
    if index == 1
        a.∂²ϵ_∂ρ²
    elseif index == 2
        a.∂²ϵ_∂ρ∂σ
    elseif index == 3
        a.∂²ϵ_∂σ²
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" Third derivative from GGA

Include the third derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
immutable GGAThirdDerivative{T0, T1, T2, T3} <: NamedTuple
    ∂³ϵ_∂ρ³::T0
    ∂³ϵ_∂ρ²∂σ::T1
    ∂³ϵ_∂ρ∂σ²::T3
    ∂³ϵ_∂σ³::T3
end
Base.length(a::GGAThirdDerivative) = 4
function Base.getindex(a::GGAThirdDerivative, index::Integer)
    if index == 1
        a.∂³ϵ_∂ρ³
    elseif index == 2
        a.∂³ϵ_∂ρ²∂σ
    elseif index == 3
        a.∂³ϵ_∂ρ∂σ²
    elseif index == 3
        a.∂³ϵ_∂σ³
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" Holds GGA energy and first derivatives """
immutable GGAEnergyAndPotential{T0, T1, T2} <: NamedTuple
    ϵ::T0
    ∂ϵ_∂ρ::T1
    ∂ϵ_∂σ::T2
end
Base.length(a::GGAEnergyAndPotential) = 3
function Base.getindex(a::GGAEnergyAndPotential, index::Integer)
    if index == 1
        a.ϵ
    elseif index == 2
        a.∂ϵ_∂ρ
    elseif index == 3
        a.∂ϵ_∂σ
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

""" All outputs from LDA """
immutable AllGGA{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9} <: NamedTuple
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
Base.length(a::AllGGA) = 10
function Base.getindex(a::AllGGA, index::Integer)
    if index == 1
        a.ϵ
    elseif index == 2
        a.∂ϵ_∂ρ
    elseif index == 3
        a.∂ϵ_∂σ
    elseif index == 4
        a.∂²ϵ_∂ρ²
    elseif index == 5
        a.∂²ϵ_∂ρ∂σ
    elseif index == 6
        a.∂²ϵ_∂σ²
    elseif index == 7
        a.∂³ϵ_∂ρ³
    elseif index == 8
        a.∂³ϵ_∂ρ²∂σ
    elseif index == 9
        a.∂³ϵ_∂ρ∂σ²
    elseif index == 10
        a.∂³ϵ_∂σ³
    else
        throw(BoundsError("Index $index is out of bounds"))
    end
end

Base.start(iter::NamedTuple) = 1
Base.next(iter::NamedTuple, state::Integer) = iter[state], state + 1
Base.done(iter::NamedTuple, state::Integer) = state > length(iter)
