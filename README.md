* A handy utility package to go along with R package CRF (https://cran.r-project.org/web/packages/CRF/index.html). 

* To install:

1. First install devtools, BiocManager, graph, Rgraphviz and RBGL:

	install.packages("devtools")
	install.packages("BiocManager")
	BiocManager::install("graph", force = TRUE)
	BiocManager::install("Rgraphviz", force = TRUE)
	BiocManager::install("RBGL", force = TRUE)

2. Turn on devtools and install CRFutil from github. CRFutil should automatically install CRF and gRbase:

	library(devtools)

	install_github("npetraco/CRFutil")

3. Test and see if the CRFutil library loads. No error messages is a sign of success. Don't worry about any warning messages.

	library(CRFutil)

* Extensive notes and tutorials on the use of this code with pair-wise Markov Random Fields at:

	https://jjcweb.jjay.cuny.edu/npetraco/tutorials/R/MRF/Notes.pptx