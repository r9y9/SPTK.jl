# extend vector to vector transformation for matrix input
for f in [:mcep,
          :gcep,
          :mgcep,
          :uels,
          :fftcep,
          :mfcc,
          :mc2b,
          :b2mc,
          :c2ir,
          :c2ndps,
          :ndps2c,
          :gc2gc,
          :lpc,
          :lpc2c,
          :lpc2lsp,
          :lpc2par,
          :lsp2sp,
          :gnorm,
          :ignorm,
          :freqt,
          :frqtr,
          :mgc2mgc,
          :mgc2sp,
          :mgclsp2sp,
          ]
    @eval begin
        function $f(x::Matrix{Float64}, args...; kargs...)
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
