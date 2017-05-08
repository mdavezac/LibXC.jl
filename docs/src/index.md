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

The functionals expect input arrays ρ and ∇ρ, and (optionally) a number of output
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

  |GGA        | unpolarized | polarized                 |
  |-----------|-------------|---------------------------|
  |ρ          | any         | `(2, ...)`                |
  |∇ρ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |ϵ          | `size(ρ)`   | `size(ρ)[2:end]`          |
  |∂ϵ/∂ρ      | `size(ρ)`   | `size(ρ)`                 |
  |∂ϵ/∂∇ρ     | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂ρ²    | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂ρ∂∇ρ  | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂²ϵ/∂∇ρ²   | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ³    | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ²∂∇ρ | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
  |∂³ϵ/∂ρ∂∇ρ² | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
  |∂³ϵ/∂∇ρ³   | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |

For the exact meaning of each dimension in each array, please refer to
[libxc](http://octopus-code.org/wiki/Libxc)

## A word about physical units

The underlying C library expects inputs in Hartree atomic units. It is possible (and
recommended) to make units part of the type of the inputs and outputs. When using Hartree
atomic units with `Cdouble`, as shown below, this will not incur any overhead. We use
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

Once a functional is created, it can be called with a number of methods to compute the
energy, the potential, as well as the second and third energy derivatives (when available
for that functional).

```jldoctest
julia> func = XCFunctional(:lda_x, false);

julia> energy(func, Cdouble[1, 2, 3]u"rho")
3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:
 -0.738559 Eₕ
 -0.930526 Eₕ
  -1.06519 Eₕ
```

Note that we create an array of `Cdouble` (with the right units, as well). The underlying C
library expects this type. Other types (and units, if not in Hartree atomic units) will
incur the cost of creating of a new array with the right type.

The following functions are available:

* [`energy`](@ref)
* [`potential`](@ref)
* [`energy_and_potential`](@ref)
* [`second_energy_derivative`](@ref)
* [`third_energy_derivative`](@ref)
* [`lda`](@ref) (all possible lda for the given functional outputs)
* [`gga`](@ref) (all possible gga outputs for the given functional)

All these functions have overloads which hide the creation of a functional from the user:

```jldoctest
julia> energy(:lda_x, [1, 2, 3]u"ρ")
3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:
 -0.738559 Eₕ
 -0.930526 Eₕ
  -1.06519 Eₕ

julia> energy(:lda_x, [1 2 3; 3 2 1]u"ρ")
3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:
 -1.23917 Eₕ
 -1.17239 Eₕ
 -1.23917 Eₕ

julia> energy(:lda_x, false, [1 2 3; 3 2 1]u"ρ")
2×3 Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},2}:
 -0.738559 Eₕ  -0.930526 Eₕ   -1.06519 Eₕ
  -1.06519 Eₕ  -0.930526 Eₕ  -0.738559 Eₕ
```

In most cases, the overhead of creating and destroying a C functional object at each call is
likely too small to matter.

The spin-polarization can be specified in the second argument (`true` for spin-polarized,
`false` for spin-polarized). If this argument is not given, then a best-guess attempt is
made: the functional is spin-polarized when ρ is at least two-dimensional and the *first*
dimension of ρ is two (`size(ρ, 1) == 2`), and the functional is unpolarized in all other
cases.

Finally, it is possible to give inputs in different units. However, this will incur the cost
of converting the array to the Hartree atomic units, both in terms of memory (an extra array
is allocated) and in terms of compute (the actual conversion). The return is always in
atomic units:

```jldoctest
julia> energy(:lda_x, false, [1 2 3; 3 2 1]u"nm^-3")
2×3 Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},2}:
 -0.0390828 Eₕ  -0.0492413 Eₕ  -0.0563672 Eₕ
 -0.0563672 Eₕ  -0.0492413 Eₕ  -0.0390828 Eₕ
```

## Using pre-allocated output array

Similar functions exist that take pre-allocated output arrays. Following `Julia`
conventions, these functions are named `energy!`, `potential!`, etc... Each function named
above has an `xxx!` counterpart.

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
(ϵ = [-0.738559,-0.930526,-1.06519]u"Eₕ", ∂ϵ_∂ρ = [-0.984745,-1.2407,-1.42025]u"∂ϵ_∂ρ")

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

