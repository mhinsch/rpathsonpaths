PACKAGE = rpathsonpaths
RSCRIPT = Rscript --no-init-file

all: install


R/RcppExports.R: src/dir_network.h
	${RSCRIPT} -e 'library(Rcpp); compileAttributes(".")'

README.md: README.Rmd
	${RSCRIPT} -e 'library(rmarkdown); render("README.Rmd", output_format="github_document")'

man/popsnetwork.Rd: R/RcppExports.R vignettes
	${RSCRIPT} -e "library(roxygen2); roxygenize()"

mans:
	${RSCRIPT} -e "library(roxygen2); roxygenize()"

vignette: 
	${RSCRIPT} -e 'library(rmarkdown); render("vignettes/overview.Rmd")'

html:
	${RSCRIPT} -e 'pkgdown::build_site()'

docs: mans vignette html

cdocs:
	doxygen

install: man/popsnetwork.Rd 
	R CMD INSTALL .

build: man/popsnetwork.Rd README.md
	R CMD build .

check: build
	R CMD check `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`

clean:
	rm -f src/*.o src/libpathsonpaths/*.o
