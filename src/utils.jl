# Utils

# used in amcep in C
function phidf!(x::Cdouble, order, α, delay::StridedVector{Cdouble})
    if length(delay) != order+1
        throw(ArgumentError("inconsistent order or delay"))
    end

    ccall((:phidf, libSPTK), Void, (Cdouble, Cint, Cdouble, Ptr{Cdouble}),
          x, order, α, delay)
end

function lspcheck(lsp::StridedVector{Cdouble})
    ccall((:lspcheck, libSPTK),
          Cint, (Ptr{Cdouble}, Cint), lsp, length(lsp)-1)
end
