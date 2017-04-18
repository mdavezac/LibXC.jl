```@meta
CurrentModule = LibXC
DocTestSetup = quote
    using LibXC
    using Unitful
    using UnitfulHartree
    func = XCFunctional(:lda_x, false)
end
```
# LibXC

This package brings julia bindings for [libxc](http://octopus-code.org/wiki/Libxc).
Currently, only bindings for LDA and GGA functionals exist. Bindings for the routines
returning information about the functional, e.g. a longer [`description`](@ref) or
[`citations`](@ref) are also given.

## Creating a Functional

The most efficient way to access the functionals is to create one:

```jldoctest
julia> functional = XCFunctional(:lda_x, true)
LibXC.XCFunctional{Float64}(:lda_x, polarized)
  - description: Slater exchange
  - kind: exchange
  - family: lda
  - spin: polarized
  - citations:
      * P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)
      * F. Bloch, Z. Phys. 57, 545 (1929)
```

The first argument is the name of the functional, and the second (optional, defaults to
true) is whether the functional is spin-polarized. The names of the available functionals
can be obtained with

```@docs
libxc_functionals
```

The functionals can be queried for their [`kind`](@ref), [`family`](@ref),
[`description`](@ref), [`citations`](@ref), and [`spin`](@ref).

## A word about input dimensionality

The functionals expect input arrays ρ and ∇ρ=|∇ρ|, and (optionally) a number of output
arrays, for the energy `ϵ`, and the different derivatives, e.g. ∂ϵ/∂ρ. Because we are
accessing a C library, some care must be taken when creating these arrays.

* All arrays must be dense (contiguous in memory) and the element type must match `Cdouble`
  (On most systems, `Cdouble` is an alias for `Float64`).
* non-spin-polarized cases: ρ can have any dimensions. All other arrays must match ρ.
* spin-polarized cases: the **first** dimension of ρ must be 2: `size(ρ) = (2, ....)`. Other
  arrays must match in very specific ways:

  |LDA       | unpolarized | polarized                 |
  |----------|-------------|---------------------------|
  |ρ         | any         | `(2, ...)`                |
  |ϵ         | `size(ρ)`   | `size(ρ)[2:end]`          |
  |∂ϵ/∂ρ     | `size(ρ)`   | `size(ρ)`                 |
  |∂²ϵ/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ³   | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |

  |GGA       | unpolarized | polarized                 |
  |----------|-------------|---------------------------|
  |ρ         | any         | `(2, ...)`                |
  |∇ρ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |ϵ         | `size(ρ)`   | `size(ρ)[2:end]`          |
  |∂ϵ/∂ρ     | `size(ρ)`   | `size(ρ)`                 |
  |∂ϵ/∂∇ρ     | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂ρ∂∇ρ  | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂∇ρ²   | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ³   | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ²∂∇ρ | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ∂∇ρ² | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
  |∂³ϵ/∂∇ρ³   | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |

For the exact meaning of each dimension in each array, please refer to
[libxc](http://octopus-code.org/wiki/Libxc)

## A word about physical units

The underlying C library expects inputs in Hartree atomic units. It is possible (and
recommended) to make units part of the type of the inputs and outputs. We use
[Unitful](http://ajkeller34.github.io/Unitful.jl/stable/),
[UnitfulHartree](https://github.com/mdavezac/UnitfulHartree.jl), and to defined (within
`LibXC.DFTUnits`) a set of units to represent the electronic density, its gradient, the
exchange-correlation energy densities, and their derivatives. These units can be accessed in
the standard way:

```jldoctest
julia> using LibXC;

julia> 1u"ρ"
1 ρ

julia> 1u"∇ρ"
1 ∇ρ

julia> 1u"grho"
1 ∇ρ

julia> 1u"ϵ"
1 ϵ

julia> 1u"Exc"
1 ϵ

julia> 1u"∂²ϵ_∂ρ∂∇ρ"
1 ∂²ϵ_∂ρ∂∇ρ
```

ρ, ∇ρ (gradient of ρ) and ϵ have non-unicode aliases, for ease of access. The energy
derivatives do not.

## Using the functionals

The majority of the functionality is contained in four functions, [`energy`](@rf),
[`potential`](@ref), [`second_energy_derivative`](@ref),
[`third_energy_derivative`](@ref). Note that user-friendly overloads exist that make
pre-allocating a functional. However, those options do come with some small overhead with
each call.

```jldoctest
julia> func = XCFunctional(:lda_x, false);

julia> energy(func, Cdouble[1, 2, 3]u"rho")
3-element Array{Quantity{Float64, Dimensions:{𝐄^-1 𝐋^2 𝐌 𝐓^-2}, Units:{ϵ}},1}:
 -0.738559 ϵ
 -0.930526 ϵ
  -1.06519 ϵ
```

Note that we create an array of `Cdouble` (with the right units, as well). The underlying C
library expects this type. Other types (and units, if not in Hartree atomic units) will
incurr the cost of creating of a new array with the right type. More complicated functions
will modify existing arrays, thus removing inefficiencies due to memory allocation:

```@meta
DocTestSetup = quote
    using LibXC
    using Unitful
    func = XCFunctional(:lda_x, false)
end
```

```jldoctest
julia> ρ = Cdouble[1, 2, 3]u"rho";

julia> ϵ = similar(ρ, LibXC.Units.ϵ{Cdouble});

julia> ∂ϵ_∂ρ = similar(ρ, LibXC.Units.∂ϵ_∂ρ{Cdouble});

julia> result = energy_and_potential!(func, ρ, ϵ, ∂ϵ_∂ρ)
(ϵ = [-0.738559,-0.930526,-1.06519]u"ϵ", ∂ϵ_∂ρ = [-0.984745,-1.2407,-1.42025]u"∂ϵ_∂ρ")

julia> result.ϵ === ϵ
true
```

For convenience, some of the functions with more complex outputs return a named tuple.
However, notice that the arrays in the tuple are aliases to the input arrays.


## API

```@autodocs
Modules = [LibXC]
```

## Available LDA functionals

```@eval
using Base.Markdown
using LibXC
result = Base.Markdown.MD()
for (name, key) in sort([u for u in LibXC.FUNCTIONALS], by=x -> x[2])
    if family(name) == LibXC.Constants.lda
        push!(result.content, LibXC.to_markdown(name))
    end
end
result
```

## Available GGA functionals

```@eval
using Base.Markdown
using LibXC
result = Base.Markdown.MD()
for (name, key) in sort([u for u in LibXC.FUNCTIONALS], by=x -> x[2])
	if family(name) == LibXC.Constants.gga
		push!(result.content, LibXC.to_markdown(name))
    end
end
result
```

