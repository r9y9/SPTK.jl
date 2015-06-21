# Utils

# used in amcep in C
function phidf!(x::Cdouble, order, α, delay::Vector{Cdouble})
    if length(delay) != order+1
        throw(ArgumentError("inconsistent order or delay"))
    end

    ccall((:phidf, libSPTK), Void, (Cdouble, Cint, Cdouble, Ptr{Cdouble}),
          x, order, α, delay)
end

function lspcheck(lsp::Vector{Cdouble})
    ccall((:lspcheck, libSPTK),
          Cint, (Ptr{Cdouble}, Cint), lsp, length(lsp)-1)
end
