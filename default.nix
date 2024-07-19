with import <nixpkgs> {};
let
  rlibs = with rPackages; [
  # base
    R

  # Imports
    # R
    data_table
    fs
    Matrix
    parallelly

    # Rcpp
    Rcpp
    RcppArmadillo
    RcppProgress

  # Suggests
    biscuiteer
    bsseq
    DelayedArray
    GenomicRanges

  # dev dependencies
    bench
    devtools
    lobstr
    pkgbuild
    pkgdown
    RcppClock
    roxygen2
    styler
    usethis
    testthat

    # RCMD check
    V8
    rhub
  ];

  _libs = with pkgs; [
  # pkg deps
    htslib
  # R CMD check
    html-tidy
    texlive.combined.scheme-full
    checkbashisms
  # dev
    ccls
    pkg-config
  ];

  in mkShell {
    nativeBuildInputs = [
      rlibs
      _libs
    ];
    shellHook = ''
      export I_R=${pkgs.R}/lib/R/include/
      export I_RCPP=${rPackages.Rcpp}/library/Rcpp/include/
      export I_ARMA=${rPackages.RcppArmadillo}/library/RcppArmadillo/include/
      export I_CLOCK=${rPackages.RcppClock}/library/RcppClock/include/
      export I_PROGRESS=${rPackages.RcppProgress}/library/RcppProgress/include/

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
      EOF
      '';
  }

# vim:set et sw=2 ts=2:
