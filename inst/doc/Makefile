RSCRIPT=Rscript

all:
	${RSCRIPT} -e "rmarkdown::render('GSE39672.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE39672.Rmd', output_format = BiocStyle::pdf_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE52980.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE52980.Rmd', output_format = BiocStyle::pdf_document2())"
	${RSCRIPT} -e "rmarkdown::render('omicsPrint.Rmd')"

