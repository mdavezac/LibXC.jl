__precompile__()
module LibXC
export description, kind, family, flags, citations, spin, energy, energy!
export potential, potential!, second_energy_derivative, third_energy_derivative
export energy_and_potential, energy_and_potential!, lda!, lda, XCFunctional, gga, gga!
export libxc_functionals, DFTUnits

using DocStringExtensions
using Unitful: Quantity, @u_str
export @u_str

macro lintpragma(s) end

if isfile(joinpath(Pkg.dir("LibXC"),"deps","deps.jl"))
    include(joinpath(Pkg.dir("LibXC"),"deps","deps.jl"))
else
    error("LibXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

""" Floating points LibXC might now about """
const CReal = Union{Cfloat, Cdouble}

include("DFTUnits.jl")
using .DFTUnits
include("units.jl")
using .Units

include("structures.jl")
include("constants.jl")
using .Constants

""" Holds citation data """
struct Citation
    ref::String
    doi::String
    bibtex::String
end
Base.show(io::IO, x::Citation) = show(io, x.ref)



const LIB_VERSION = let
    version = Cint[0, 0, 0]
    @lintpragma("Ignore use of undeclared variable libxc")
    @lintpragma("Ignore use of undeclared variable ccall")
    @eval ccall((:xc_version, $libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}),
                Ref($version), Ref($version, 2), Ref($version, 3))
    VersionNumber(version...)
end



""" Functional names and keys """
const FUNCTIONALS = let
    tags = @eval cglobal(
        (:xc_functional_keys, $(LibXC.libxc)), LibXC.CFunctionalKey)
    result = Dict{Symbol, Int32}()

    i, key = 1, unsafe_load(tags)
    while key.key > 0
        n = findfirst(x -> x == UInt32('\0'), key.name) - 1
        push!(result, Symbol(join(Char(key.name[u]) for u in 1:n)) => key.key)
        i += 1
        key = unsafe_load(tags, i)
    end
    result
end
"""" Functional keys and names """
const iFUNCTIONALS = Dict(v => k for (k, v) in FUNCTIONALS)

""" Names of the available libxc functionals """
function libxc_functionals()
    collect(filter(keys(FUNCTIONALS)) do x
        local name = string(x)
        length(name) > 3 && (name[1:3] == "lda" || name[1:3] == "gga")
    end)
end

abstract type AbstractXCFunctional end
abstract type AbstractLibXCFunctional{T <: CReal} <: AbstractXCFunctional end

""" Manages pointer to C libxc funtional data """
type XCFunctional{T <: CReal} <: AbstractLibXCFunctional{T}
    c_ptr::Ptr{CFuncType{T}}
    XCFunctional{T}(ptr::Ptr{CFuncType{T}}) where T <: CReal = new(ptr)
end

""" Creates a functional from it's name and polarization

The name should be a one of the following symbols:

- LDA: $(join((u for u in keys(FUNCTIONALS) if string(u)[1:3] == "lda"), ", "))
- GGA: $(join((u for u in keys(FUNCTIONALS) if string(u)[1:3] == "gga"), ", "))

The second argument is `true` if the functional should be polarized.
"""
function XCFunctional(name::Symbol, spin_polarized::Bool=true)
    name ∉ keys(FUNCTIONALS) && error("Functional $name does not exist")
    ptr = ccall((:xc_func_alloc, libxc), Ptr{CFuncType{Cdouble}}, ())
    functional = XCFunctional{Cdouble}(ptr)
    err = ccall((:xc_func_init, libxc), Cint, (Ptr{CFuncType{Cdouble}}, Cint, Cint),
                ptr, FUNCTIONALS[name], spin_polarized ? 2: 1)
    err ≠ 0 && error("Error $err encountered in LibXC")
    finalizer(functional, _delete_libxc_functional)
    functional
end
function XCFunctional(name::Symbol, spin::Constants.SPIN)
    XCFunctional(name, spin == Constants.polarized)
end
XCFunctional(n::Symbol, spin::Integer) = XCFunctional(n, convert(Constants.SPIN, spin))

function _func_info(func::AbstractLibXCFunctional{Cdouble})
    func_type = ccall((:xc_func_get_info, libxc), Ptr{CFuncInfoType{Cdouble}},
                      (Ptr{CFuncType{Cdouble}}, ), func.c_ptr)
    unsafe_load(func_type)
end

function _delete_libxc_functional(func::AbstractLibXCFunctional{Cdouble})
    if func.c_ptr ≠ C_NULL
        ccall((:xc_func_end, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
        ccall((:xc_func_free, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
    end
end


""" Long name/description of a functional """
description(info::CFuncInfoType) = unsafe_string(info.name)
""" Integer key of the functional """
libkey(info::CFuncInfoType) = info.number
""" Whether Exchange, Correlation, or Exchange-Correlation """
kind(info::CFuncInfoType) = convert(Constants.KIND, info.kind)
""" Whether LDA, GGA, etc """
family(info::CFuncInfoType) = convert(Constants.FAMILY, info.family)
""" Spin polarization """
spin(func::CFuncType) = convert(Constants.SPIN, func.nspin)

""" Set of flags describing the functional """
function flags(func::CFuncInfoType)
    result = Set{Constants.FLAGS}()
    for flag in Base.instances(Constants.FLAGS)
        (convert(Cint, flag) & func.flags) ≠ 0 && push!(result, flag)
    end
    result
end

""" List of journal references """
function citations(info::CFuncInfoType)
    result = Citation[]
    for ptr in info.refs
        ptr == C_NULL && break
        ref = unsafe_load(ptr)
        push!(result, Citation(unsafe_string(ref.ref), unsafe_string(ref.doi),
                               unsafe_string(ref.bibtex)))
    end
    result
end

description(func::AbstractLibXCFunctional) = description(_func_info(func))
libkey(func::AbstractLibXCFunctional) = libkey(_func_info(func))
kind(func::AbstractLibXCFunctional) = kind(_func_info(func))
family(func::AbstractLibXCFunctional) = family(_func_info(func))
flags(func::AbstractLibXCFunctional) = flags(_func_info(func))
citations(func::AbstractLibXCFunctional) = citations(_func_info(func))
spin(func::AbstractLibXCFunctional) = convert(Constants.SPIN, unsafe_load(func.c_ptr).nspin)
citations(name::Symbol) = citations(XCFunctional(name))
description(name::Symbol) = description(XCFunctional(name))
libkey(name::Symbol) = libkey(XCFunctional(name))
kind(name::Symbol) = kind(XCFunctional(name))
family(name::Symbol) = family(XCFunctional(name))
flags(name::Symbol) = flags(XCFunctional(name))

"""
    $(SIGNATURES)

Helps determine size of an output. `dims` refers to the size of input density array, and
`factor` is the size we expect for the first dimension. It will vary with the functional
(LDA, GGA, ...) and the kind of output (energy, potential, ...).
"""
function output_size(polarized::Bool, dims::NTuple, factor::Integer)
    if length(dims) == 0
        throw(ArgumentError("Empty size tuple"))
    elseif !polarized
        dims
    elseif length(dims) == 1 && polarized
        if dims[1] % 2 ≠ 0
            throw(ArgumentError("Odd array size for polarized functional"))
        end
        warn("Spin polarized function, but dimensionality of ρ is 1")
        (factor * dims[1] / 2,)
    elseif dims[1] == 2
        if factor == 1
            dims[2:end]
        else
            (factor, dims[2:end]...)
        end
    else
        throw(ArgumentError("Spin polarization expects size(ρ, 1) == 2"))
    end
end
function output_size(func::AbstractLibXCFunctional, ρ::DenseArray, factor::Integer)
    output_size(spin(func), Base.size(ρ), factor)
end
function output_size(s::Constants.SPIN, dims::NTuple, factor::Integer)
    output_size(s == Constants.polarized, dims, factor)
end


include("Checks.jl")
using .Checks: *

include("named_tuples.jl")
include("lda.jl")
include("gga.jl")
include("overloads.jl")


""" Prints functional to markdown, mostly for docs """
function to_markdown(name)
    func = XCFunctional(name, false)
    result = Base.Markdown.MD()
    push!(result.content, Base.Markdown.Header{3}(description(func)))

    list = Base.Markdown.List()
    push!(list.items,
		Any[Base.Markdown.Paragraph(["name: ", Base.Markdown.Code("", ":$name")])])
    push!(list.items, Any[Base.Markdown.Paragraph(["kind: $(kind(func))"])])

    availability = String[]
    Constants.exc ∈ flags(func) && push!(availability, "εxc")
    Constants.vxc ∈ flags(func) && push!(availability, "∂εxc")
    Constants.fxc ∈ flags(func) && push!(availability, "∂²εxc")
    Constants.fxc ∈ flags(func) && push!(availability, "∂³εxc")
    if length(availability) > 0
        push!(list.items,
              Any[Base.Markdown.Paragraph(
                    ["available derivatives: ", join(availability, ", ")])])
    end

    references = Base.Markdown.List()
    for citation in citations(func)
        link = Base.Markdown.Link(citation.ref, "https://dx.doi.org/$(citation.doi)")
        push!(references.items, Any[Base.Markdown.Paragraph([link])])
    end
    push!(list.items, Any[Base.Markdown.Paragraph(["references:"]), references])
    push!(result.content, list)
    result
end

function Base.showcompact(io::IO, func::AbstractLibXCFunctional)
    symb = iFUNCTIONALS[libkey(func)]
    print(io, "$(typeof(func))(:$symb, $(LibXC.spin(func)))")
end
function Base.show(io::IO, func::AbstractLibXCFunctional)
    showcompact(io, func)
    println(io,
        "\n  - description: $(description(func))",
        "\n  - kind: $(kind(func))",
        "\n  - family: $(family(func))",
        "\n  - spin: $(spin(func))",
        "\n  - citations:")
    for cit in citations(func)
        println(io, "      * $(cit.ref)")
    end
end



end # module
