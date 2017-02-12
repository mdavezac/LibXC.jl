module LibXC

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

function _get_version()
    version = Cint[0, 0, 0]
    ccall((:xc_version, libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}), Ref(version),
          Ref(version, 2), Ref(version, 3))
    VersionNumber(version...)
end

const lib_version = _get_version()

_xc_key(name::AbstractString) = ccall((:xc_functional_get_number), Cint, (Cstring,), name)
function _xc_name(key::Integer)
    name = ccall((:xc_functional_get_name), Cstring, (Cint,), key)
    unsafe_string(name)
end

end # module
