module Internals
export description, kind, family, flags, citations, spin, libxc_functionals,
       LIB_VERSION, FUNCTIONALS, XCFunctional, AbstractLibXCFunctional, AbstractXCFunctional
using ..Constants
using DFTShims: SpinDegenerate, ColinearSpin

macro lintpragma(s) end

if isfile(joinpath(Pkg.dir("LibXC"),"deps","deps.jl"))
    include(joinpath(Pkg.dir("LibXC"),"deps","deps.jl"))
else
    error("LibXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

""" Floating points LibXC might now about """
const CReal = Union{Cfloat, Cdouble}

include("structures.jl")

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
    tags = @eval cglobal((:xc_functional_keys, $(libxc)), CFunctionalKey)
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
libxc_functionals() = 
    collect(filter(keys(FUNCTIONALS)) do x
        local name = string(x)
        length(name) > 3 && (name[1:3] == "lda" || name[1:3] == "gga")
    end)

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
XCFunctional(name::Symbol, spin_polarized::Bool=true) = begin
    name ∉ keys(FUNCTIONALS) && error("Functional $name does not exist")
    ptr = ccall((:xc_func_alloc, libxc), Ptr{CFuncType{Cdouble}}, ())
    functional = XCFunctional{Cdouble}(ptr)
    err = ccall((:xc_func_init, libxc), Cint, (Ptr{CFuncType{Cdouble}}, Cint, Cint),
                ptr, FUNCTIONALS[name], spin_polarized ? 2: 1)
    err ≠ 0 && error("Error $err encountered in LibXC")
    finalizer(functional, _delete_libxc_functional)
    functional
end
XCFunctional(name::Symbol, spin::Constants.SPIN) =
    XCFunctional(name, spin == Constants.polarized)
XCFunctional(n::Symbol, spin::Integer) = XCFunctional(n, convert(Constants.SPIN, spin))
XCFunctional(n::Symbol, ::SpinDegenerate) = XCFunctional(n, false)
XCFunctional(n::Symbol, ::ColinearSpin) = XCFunctional(n, true)

_func_info(func::AbstractLibXCFunctional{Cdouble}) = begin
    func_type = ccall((:xc_func_get_info, libxc), Ptr{CFuncInfoType{Cdouble}},
                      (Ptr{CFuncType{Cdouble}}, ), func.c_ptr)
    unsafe_load(func_type)
end

_delete_libxc_functional(func::AbstractLibXCFunctional{Cdouble}) = begin
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
flags(func::CFuncInfoType) = begin
    result = Set{Constants.FLAGS}()
    for flag in Base.instances(Constants.FLAGS)
        (convert(Cint, flag) & func.flags) ≠ 0 && push!(result, flag)
    end
    result
end

""" List of journal references """
citations(info::CFuncInfoType) = begin
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

""" Prints functional to markdown, mostly for docs """
to_markdown(name) = begin
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

Base.showcompact(io::IO, func::AbstractLibXCFunctional) = begin
    symb = iFUNCTIONALS[libkey(func)]
    print(io, "$(typeof(func))(:$symb, $(spin(func))")
end
Base.show(io::IO, func::AbstractLibXCFunctional) = begin
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
end
