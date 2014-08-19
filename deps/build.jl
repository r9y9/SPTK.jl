using BinDeps

@BinDeps.setup

deps = [
        sptk = library_dependency("libSPTK")
        ]
        
const sptkversion = "3.7.0"

provides(Sources,
         URI("https://github.com/r9y9/SPTK/archive/v$(sptkversion).tar.gz"),
         sptk,
         unpacked_dir="SPTK-$(sptkversion)")
         
prefix=joinpath(BinDeps.depsdir(sptk), "usr")
sptksrcdir = joinpath(BinDeps.depsdir(sptk),"src", "SPTK-$(sptkversion)")

provides(SimpleBuild,
          (@build_steps begin
              GetSources(sptk)
              @build_steps begin
                  ChangeDirectory(sptksrcdir)
                  `./waf configure --prefix=$prefix`
                  `./waf build`
                  `./waf install`
              end
           end), sptk, os = :Unix)

@BinDeps.install [:libSPTK => :libSPTK]
