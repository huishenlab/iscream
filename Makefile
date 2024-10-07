all: clean build
docs: roxygen clean build

build:
	R -e "Rcpp::compileAttributes()"
	R CMD INSTALL .

roxygen:
	R -e "Rcpp::compileAttributes(); roxygen2::roxygenise()"

clean:
	rm -rf src/*.o src/*.so iscream_0.0.0.9000.tar.gz iscream.Rcheck

check:
	R CMD build . && R CMD check *.tar.gz

cclean:
	rm -rf iscream_0.0.0.9000.tar.gz iscream.Rcheck

site:
	R -e "pkgdown::build_site()"

.PHONY: all build clean
