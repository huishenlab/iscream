{
  description = "Flake to get iscream development environment";
  inputs.nixpkgs.url = "github:rstats-on-nix/nixpkgs/r-bioc-devel";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};

      LinkingTo = with pkgs.rPackages; [
        Rcpp
        RcppArmadillo
        RcppProgress
        RcppSpdlog
        Rhtslib
        stringfish
      ];

      Imports = with pkgs.rPackages; [
        data_table
        Matrix
        parallelly
      ];

      Suggests = with pkgs.rPackages; [
        BiocStyle
        GenomicRanges
        ggplot2
        ggridges
        SummarizedExperiment
      ];

      rDevDeps = with pkgs.rPackages; [
        BiocCheck
        BiocManager
        BiocVersion
        covr
        devtools
        DT
        htmltools
        lobstr
        pkgdown
        rhub
        roxygen2
        styler
        testthat
        usethis
        V8
      ];

      Bioc = with pkgs.rPackages; [
        biscuiteer
        bsseq
      ];

      htslib = pkgs.htslib.overrideAttrs (finalAttrs: previousAttrs: {
        buildInputs = previousAttrs.buildInputs ++ [ pkgs.libdeflate ];
      });
      sysDeps = with pkgs; [
        R
        gcc
        htslib
        pkg-config
      ];

      sysDevDeps = with pkgs; [
        air-formatter
        html-tidy
        texlive.combined.scheme-full
        checkbashisms
        ccls
      ];

      # default package
      rDeps = [ LinkingTo Imports Suggests ];
      iscream = pkgs.rPackages.buildRPackage {
        name = "iscream";
        src = self;
        nativeBuildInputs = sysDeps;
        propagatedBuildInputs = rDeps;
      };
      # Create R development environment with iscream and other useful libraries
      rvenv = pkgs.rWrapper.override {
        packages = rDeps ++ rDevDeps ++ sysDeps ++ sysDevDeps;
      };
    in {
      packages.default = iscream;
      devShells.default = pkgs.mkShell {
          buildInputs = rDeps ++ rDevDeps ++ sysDeps ++ sysDevDeps;
          inputsFrom = pkgs.lib.singleton iscream;
          packages = pkgs.lib.singleton rvenv;
          shellHook = ''
            export I_R=${pkgs.R}/lib/R/include/
            export I_RCPP=${pkgs.rPackages.Rcpp}/library/Rcpp/include/
            export I_ARMA=${pkgs.rPackages.RcppArmadillo}/library/RcppArmadillo/include/
            export I_PROGRESS=${pkgs.rPackages.RcppProgress}/library/RcppProgress/include/
            export I_LOG=${pkgs.rPackages.RcppSpdlog}/library/RcppSpdlog/include/
            export I_STRINGFISH=${pkgs.rPackages.stringfish}/library/stringfish/include/

            export I_HTSLIB=${pkgs.htslib}/include/
            export L_HTSLIB=${pkgs.htslib}/lib/libhts.a
            export L_CURL=${pkgs.curl.out}/lib/libcurl.so

            mkdir -p "$HOME/.R"
            export R_LIBS_USER="$HOME/.R"

            cat > .ccls << EOF
            clang
            %c -std=c11
            %cpp -std=c++2a
            -I$I_R
            -I$I_RCPP
            -I$I_ARMA
            -I$I_HTSLIB
            -I$I_CLOCK
            -I$I_PROGRESS
            -I$I_LOG
            -I$I_STRINGFISH
            EOF
            '';
      };
    });
}
