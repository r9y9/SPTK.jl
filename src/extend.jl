# List of vector to vector conversion functions
# should accept a function call like: f(x) where x is a input vector
const vec2vec = [
                 :mcep,
                 :gcep,
                 :mgcep,
                 :uels,
                 :fftcep,
                 :lpc,

                 :mfcc,

                 :lpc2c,
                 :lpc2lsp,
                 :lpc2par,
                 :par2lpc,
                 :lsp2sp,

                 :mc2b,
                 :b2mc,
                 :b2c,
                 :c2acr,
                 :c2ir,
                 :ic2ir,
                 :c2ndps,
                 :ndps2c,
                 :gc2gc,
                 :gnorm,
                 :ignorm,
                 :freqt,
                 :frqtr,
                 :mgc2mgc,
                 :mgc2sp,
                 :mgclsp2sp
                 ]

# extend vector to vector conversion for matrix input (col-wise)
# should accept a function call like: f(x) where x is a input matrix
for f in vec2vec
    @eval begin
        function $f(x::StridedMatrix{Cdouble}, args...; kargs...)
            outbuf = $f(view(x, :, 1), args...; kargs...)
            ret = Array{eltype(outbuf)}(length(outbuf), size(x, 2))
            copy!(ret, 1, outbuf, 1, length(outbuf))
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(view(x, :, i), args...; kargs...)
            end
            ret
        end
    end
end
