RSCRIPT=Rscript

all:
	${RSCRIPT} -e "rmarkdown::render('GSE39672.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE52980.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE102177_siblings.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('GSE100227_sisters.Rmd', output_format = BiocStyle::html_document2())"
	${RSCRIPT} -e "rmarkdown::render('omicsPrint.Rmd', output_format = BiocStyle::html_document2())"


