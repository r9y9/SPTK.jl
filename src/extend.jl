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
        function $f(x::Matrix{Cdouble}, args...; kargs...)
            # TODO: avoid allocations (pass SubArray instead of Array?)
            r = $f(x[:, 1], args...; kargs...)
            ret = Array(eltype(r), size(r, 1), size(x, 2))
            for i = 1:length(r)
                @inbounds ret[i, 1] = r[i]
            end
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(x[:, i], args...; kargs...)
            end
            ret
        end
    end
end
