# Supplementary Material
The code provided by this repository is part of the supplementary material for the scientific article "Comparison of Three Untargeted Data Processing Workflows for Evaluating LC-HRMS Metabolomics Data" by Selina Hemmer et al., Metabolites, 2020 (DOI:10.3390/metabo10090378). See the following sections for detailed information about it.

## Data Files
The folder "mzXMl\_data\_files" contains all the analyses that were used to evaluation the workflows. There are two subfolders containing measurements for each ionization mode.

## Source Code
The source code can be found within the folder "src". Refer to the following section to learn about how to use it.

### Dependencies
The code deployed by this repository requires the following packages:
- tidyverse (CRAN)
- ggrepel (CRAN)
- gplots (CRAN)
- MASS (CRAN)
- Rtsne (CRAN)
- XCMS (Bioconductor)
- CAMERA (Bioconductor)
- metaboLib.R (Github)

If not installed you can use the following lines to install them:

	install.packages(c("tidyverse", "ggrepel", "gplots", "MASS", "Rtsne")
	source("https://bioconductor.org/biocLite.R")
	if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("xcms", "CAMERA")
	
metabolib.R can be found in its corresponding [GitHub](https://github.com/saskema/metaboLib) repository and merely needs to be placed in your working directory.

### Run Code
To run the code, place the corresponding R script into your working directory along with metaboLib.R and the mzXMl files. "HLM\_neg.R" was used for processing files obtained after analysis using negative ionization mode, "HLM\_pos.R" was used for correspoding analyses in positive ionization mode.

