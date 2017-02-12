module LibXC

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

function _get_version()
    version = Cint[0, 0, 0]
    ccall((:xc_version, libxc), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), Ref(version),
          Ref(version, 2), Ref(version, 3))
    VersionNumber(version...)
end

const version = _get_version()

_xc_key(name::AbstractString) = ccall((:xc_functional_get_number), Cint, (Cstring,), name)
_xc_name(key::Integer) = ccall((:xc_functional_get_number), Cstring, (Cint,), key)


end # module
