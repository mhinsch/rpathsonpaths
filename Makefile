PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: install


RCPP:
	${RSCRIPT} -e 'library(Rcpp); compileAttributes(".")'

README.md: README.Rmd
	${RSCRIPT} -e 'library(rmarkdown); render("README.Rmd", output_format="github_document")'

roxygen: RCPP vignettes
	@mkdir -p man
	${RSCRIPT} -e "library(roxygen2); roxygenize()"

vignettes: 
	${RSCRIPT} -e 'library(rmarkdown); render("vignettes/overview.Rmd")'

install: roxygen README.md
	R CMD INSTALL .

build: roxygen README.md
	R CMD build .

check: build
	R CMD check `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck


