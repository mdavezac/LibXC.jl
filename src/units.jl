using Unitful
using UnitfulHartree

typealias CReal Union{Cfloat, Cdouble}
typealias Density{T <: CReal} Quantity{T, typeof(dimension(u"𝖊*a₀^-3")), typeof(u"𝖊*a₀^-3")}
typealias EnergyDensity{T <: CReal} Quantity{T, typeof(dimension(u"Eₕ/𝖊")), typeof(u"Eₕ/𝖊")}
typealias EnergyDensityDerivatives{Nrho, Ngrho, T <: CReal}
    typeof(EnergyDensity(1)/Density(1)^(Nrho+Ngrho)*(1u"a₀")^Ngrho)

typealias DensityArray{T <: CReal} DenseArray{Density, N}
typealias EnergyDensityArray{T <: CReal, N} DenseArray{EnergyDensity, N}
typealias PotentialDensityArray{T <: CReal, N} DenseArray{PotentialDensity, N}

Base.one{T}(::Type{Density{T}}) = Density(one(T))
