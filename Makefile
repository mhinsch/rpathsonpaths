PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: install


RCPP:
	${RSCRIPT} -e 'library(Rcpp); compileAttributes(".")'

README:
	${RSCRIPT} -e 'library(rmarkdown); render("README.Rmd", output_format="github_document")'

roxygen: RCPP
	@mkdir -p man
	${RSCRIPT} -e "library(roxygen2); roxygenize()"

install: roxygen README
	R CMD INSTALL .

build: roxygen README
	R CMD build .

check: build
	R CMD check `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck


