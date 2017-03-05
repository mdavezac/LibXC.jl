module LibXC
export description, kind, family, flags, citations, spin, energy, energy!
export potential, potential!, second_energy_derivative, third_energy_derivative
export energy_and_potential, energy_and_potential!, lda!, lda, XCFunctional

using NamedTuples: @NT

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

include("structures.jl")
include("constants.jl")

""" Holds citation data """
immutable Citation
    ref::String
    doi::String
    bibtex::String
end
Base.show(io::IO, x::Citation) = show(io, x.ref)



const LIB_VERSION = let
    version = Cint[0, 0, 0]
    @eval ccall((:xc_version, $libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}),
                Ref($version), Ref($version, 2), Ref($version, 3))
    VersionNumber(version...)
end



""" Functional names and keys """
const FUNCTIONALS = let
    keys = @eval cglobal(
        (:xc_functional_keys, $(LibXC.libxc)), LibXC.CFunctionalKey)
    result = Dict{Symbol, Int32}()

    i, key = 1, unsafe_load(keys)
    while key.key > 0
        n = findfirst(x -> x == UInt32('\0'), key.name) - 1
        push!(result, Symbol(join(Char(key.name[u]) for u in 1:n)) => key.key)
        i += 1
        key = unsafe_load(keys, i)
    end
    result
end

abstract AbstractXCFunctional
abstract AbstractLibXCFunctional{T <: CReal} <: AbstractXCFunctional

type XCFunctional{T <: CReal} <: AbstractLibXCFunctional{T}
    c_ptr::Ptr{CFuncType{T}}
    XCFunctional(ptr::Ptr{CFuncType{T}}) = new(ptr)
end

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

function Base.size(polarized::Bool, dims::NTuple, factor::Integer)
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
function Base.size(func::AbstractLibXCFunctional, ρ::DenseArray, factor::Integer)
    Base.size(spin(func), size(ρ), factor)
end
function Base.size(s::Constants.SPIN, dims::NTuple, factor::Integer)
    Base.size(s == Constants.polarized, dims, factor)
end
function Base.size(s::Union{Bool, Constants.SPIN}, ρ::DenseArray, factor::Integer)
    Base.size(s, size(ρ), factor)
end

for (name, factor) ∈ [(:esize, 1), (:vsize, 2), (:fsize, 3), (:ksize, 4)]
    @eval begin
        function $name(polarized::Bool, dims::NTuple)
            if length(dims) == 0
                throw(ArgumentError("Empty size tuple"))
            elseif !polarized
                dims
            elseif length(dims) == 1 && polarized
                if dims[1] % 2 ≠ 0
                    throw(ArgumentError("Odd array size for polarized functional"))
                end
                warn("Spin polarized function, but dimensionality of ρ is 1")
                $(factor == 1 ? :(dims[1]/2,): :(($factor * dims[1]/2,)))
            elseif dims[1] == 2
                $(factor == 1 ? :(dims[2:end]): :(($factor, dims[2:end]...)))
            else
                throw(ArgumentError("Spin polarization expects size(ρ, 1) == 2"))
            end
        end
        $name(func::AbstractLibXCFunctional, ρ::DenseArray) = $name(spin(func), size(ρ))
        $name(s::Constants.SPIN, dims::NTuple) = $name(s == Constants.polarized, dims)
        $name(s::Union{Bool, Constants.SPIN}, ρ::DenseArray) = $name(s, size(ρ))
    end
end

include("lda.jl")
include("gga.jl")
end # module
