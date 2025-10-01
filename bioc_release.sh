#!/usr/bin/env sh

git rm -r --cached \
    vignettes/*.nocompute \
    pkgdown \
    .github \
    Makefile \
    .lintr \
    air.toml \
    flake.* \
    codecov.yml \
    bioc_release.sh

git commit -a -m "chore(bioc): create release package"
