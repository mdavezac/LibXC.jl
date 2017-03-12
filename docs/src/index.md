```@meta
CurrentModule = LibXC
DocTestSetup = quote
    using LibXC
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

The functionals expect input arrays ρ and σ=|∇ρ|, and (optionally) a number of output
arrays, for the energy `εxc`, and the different derivatives, e.g. ∂εxc/∂ρ. Because we are
accessing a C library, some care must be taken when creating these arrays.

* All arrays must be dense (contiguous in memory) and the element type must match `Cdouble`
  (On most systems, `Cdouble` is an alias for `Float64`).
* non-spin-polarized cases: ρ can have any dimensions. All other arrays must match ρ.
* spin-polarized cases: the **first** dimension of ρ must be 2: `size(ρ) = (2, ....)`. Other
  arrays must match in very specific ways:

  |LDA         | unpolarized | polarized                 |
  |------------|-------------|---------------------------|
  |ρ           | any         | `(2, ...)`                |
  |εxc         | `size(ρ)`   | `size(ρ)[2:end]`          |
  |∂εxc/∂ρ     | `size(ρ)`   | `size(ρ)`                 |
  |∂²εxc/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂³εxc/∂ρ³   | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |

  |GGA         | unpolarized | polarized                 |
  |------------|-------------|---------------------------|
  |ρ           | any         | `(2, ...)`                |
  |σ           | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |εxc         | `size(ρ)`   | `size(ρ)[2:end]`          |
  |∂εxc/∂ρ     | `size(ρ)`   | `size(ρ)`                 |
  |∂εxc/∂σ     | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²εxc/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
  |∂²εxc/∂ρ∂σ  | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂²εxc/∂σ²   | `size(ρ)`   | `(6, size(ρ)[2:end]...)`  |
  |∂³εxc/∂ρ³   | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
  |∂³εxc/∂ρ²∂σ | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
  |∂³εxc/∂ρ∂σ² | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
  |∂³εxc/∂σ³   | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |

For the exact meaning of each dimension in each array, please refer to
[libxc](http://octopus-code.org/wiki/Libxc)

## Using the functionals

The majority of the functionality is contained in four functions, [`energy`](@rf),
[`potential`](@ref), [`second_energy_derivative`](@ref),
[`third_energy_derivative`](@ref). Note that user-friendly overloads exist that make
pre-allocating a functional. However, those options do come with some small overhead with
each call.

```jldocset
julia> func = XCFunctional(:lda_x, false);

julia> func(Cdouble[1, 2, 3])
3-element Array{Float64,1}:
 -0.738559
 -0.930526
 -1.06519
```


More complicated functions will modify existing arrays, thus removing inefficiencies due to
memory allocation:

```@meta
DocTestSetup = quote
    using LibXC
    func = XCFunctional(:lda_x, false)
end
```

```jldocset
julia> ρ = Cdouble[1, 2, 3];

julia> εxc, pot = similar(ρ), similar(ρ);

julia> result = energy_and_potential(func, ρ, εxc, pot)
(energy => [-0.738559,-0.930526,-1.06519], potential => [-0.984745,-1.2407,-1.42025])

julia> result.energy === εxc
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

