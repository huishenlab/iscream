with import <nixpkgs> {};
let
  rlibs = with rPackages; [
  # base
    R

  # Imports
    fs
    Rcpp
    RcppArmadillo
    Matrix
    data_table

  # Suggests
    bsseq
    DelayedArray
    GenomicRanges

  # dev dependencies
    pkgdown
    pkgbuild
    roxygen2
  ];

  _libs = with pkgs; [
    ccls # Cpp LSP
    htslib
  ];

  in mkShell {
    nativeBuildInputs = [
      rlibs
      _libs
    ];
    shellHook = ''
      export I_RCPP=${rPackages.Rcpp}/library/Rcpp/include/
      export I_ARMA=${rPackages.RcppArmadillo}/library/RcppArmadillo/include/
      export I_HTSLIB=${pkgs.htslib}/include/
      export L_HTSLIB=${pkgs.htslib}/lib/libhts.a
      export I_R=${pkgs.R}/lib/R/include/
      export L_CURL=${pkgs.curl.out}/lib/libcurl.so

      mkdir -p "$HOME/.R"
      export R_LIBS_USER="$HOME/.R"

      cat > .ccls << EOF
      clang
      %c -std=c11
      %cpp -std=c++2a
      -I$I_RCPP
      -I$I_ARMA
      -I$I_R
      -I$I_HTSLIB
      EOF
      '';
  }

# vim:set et sw=2 ts=2:
