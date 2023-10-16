with import <nixpkgs> {};
let
  rlibs = with rPackages; [
  # base
    R

  # Imports
    fs
    Rcpp
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
  ];

  in mkShell {
    nativeBuildInputs = [
      rlibs
      _libs
    ];
    shellHook = ''
      mkdir -p "$HOME/.R"
      export R_LIBS_USER="$HOME/.R"
      '';
  }

# vim:set et sw=2 ts=2:
