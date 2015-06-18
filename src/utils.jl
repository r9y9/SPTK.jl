# Utils

# used in amcep in C
function phidf!(x::Cdouble, order, α, delay::Vector{Cdouble})
    if length(delay) != order+1
        throw(ArgumentError("inconsistent order or delay"))
    end

    ccall((:phidf, libSPTK), Void, (Cdouble, Cint, Cdouble, Ptr{Cdouble}),
          x, order, α, delay)
end

function lspcheck(lpc::Vector{Cdouble})
    r = ccall((:lspcheck, libSPTK),
              Cint, (Ptr{Cdouble}, Cint), lpc, length(lpc)-1)::Cint
    ifelse(r == 0, true, false)
end
