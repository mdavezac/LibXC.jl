module LibXC
export description, kind, family, flags, citations, spin

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
abstract AbstractLibXCFunctional <: AbstractXCFunctional

type XCFunctional <: AbstractLibXCFunctional
    c_ptr::Ptr{CFuncType{Cdouble}}
    XCFunctional(ptr::Ptr{CFuncType{Cdouble}}) = new(ptr)
end

function XCFunctional(name::Symbol, spin_polarized::Bool=true)
    name ∉ keys(FUNCTIONALS) && error("Functional $name does not exist")
    ptr = ccall((:xc_func_alloc, libxc), Ptr{CFuncType{Cdouble}}, ())
    functional = XCFunctional(ptr)
    ccall((:xc_func_init, libxc), Cint, (Ptr{CFuncType{Cdouble}}, Cint, Cint),
          ptr, FUNCTIONALS[name], spin_polarized ? 2: 1)
    finalizer(functional, _delete_libxc_functional)
    functional
end
function XCFunctional(name::Symbol, spin::Constants.SPIN)
    XCFunctional(name, spin == Constants.polarized)
end
XCFunctional(n::Symbol, spin::Integer) = XCFunctional(n, convert(Constants.SPIN, spin))

function _func_info(func::AbstractLibXCFunctional)
    func_type = ccall((:xc_func_get_info, libxc), Ptr{CFuncInfoType{Cdouble}},
                      (Ptr{CFuncType{Cdouble}}, ), func.c_ptr)
    unsafe_load(func_type)
end

function _delete_libxc_functional(func::AbstractLibXCFunctional)
    if func.c_ptr ≠ C_NULL
        ccall((:xc_func_end, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
        ccall((:xc_func_free, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
    end
end

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


end # module
