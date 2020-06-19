library(rmarkdown)

all_rmd <- list.files(path = ".", recursive = TRUE, pattern=".Rmd", full.names = TRUE)
lapply(as.list(all_rmd), render, output_dir = ".", output_format = c("html_document", "pdf_document"))