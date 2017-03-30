using Unitful
using UnitfulHartree

typealias CReal Union{Cfloat, Cdouble}
typealias Density{T <: CReal} Quantity{T, typeof(dimension(u"𝖊*a₀^-3")), typeof(u"𝖊*a₀^-3")}
typealias EnergyDensity{T <: CReal} Quantity{T, typeof(dimension(u"Eₕ/𝖊")), typeof(u"Eₕ/𝖊")}
typealias PotentialDensity{T <: CReal}
    Quantity{T, typeof(dimension(u"Eₕ*a₀^3/𝖊^2")), typeof(u"Eₕ*a₀^3/𝖊^2")}

typealias DensityArray{T <: CReal} DenseArray{Density, N}
typealias EnergyDensityArray{T <: CReal, N} DenseArray{EnergyDensity, N}
typealias PotentialDensityArray{T <: CReal, N} DenseArray{PotentialDensity, N}
