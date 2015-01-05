using BinDeps

@BinDeps.setup

deps = [
        sptk = library_dependency("libSPTK")
        ]

const version = "3.8.2"

provides(Sources,
         URI("https://github.com/r9y9/SPTK/archive/v$(version).tar.gz"),
         sptk,
         unpacked_dir="SPTK-$(version)")

prefix = joinpath(BinDeps.depsdir(sptk), "usr")
srcdir = joinpath(BinDeps.depsdir(sptk), "src", "SPTK-$(version)")

provides(SimpleBuild,
          (@build_steps begin
              GetSources(sptk)
              @build_steps begin
                  ChangeDirectory(srcdir)
                  `./waf configure --prefix=$prefix`
                  `./waf build`
                  `./waf install`
              end
           end), sptk, os = :Unix)

@BinDeps.install [:libSPTK => :libSPTK]
