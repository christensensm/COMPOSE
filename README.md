# Crispr-OMics for PhenOtypic ScrEens (COMPOSE)

Crispr-OMics for PhenOtypic ScrEens (COMPOSE) is a method developed for analyzing CRISPR screens. Functions range is usage from visualization of data for quality control to calculation of fold-changes and subsequent pathway analyses.

### Installation

COMPOSE was created using R v3.5.1 and some of the installations may not follow the same coding. Recent updates altered the way Bioconductor packages were installed for R.
Install dependencies and the package via the following:

```
install.packages("BiocManager")
BiocManager::install(c('BiocParallel','AnnotationDbi','clusterProfiler','DESeq2','fgsea','gage','pathview','ReactomePA','limma','edgeR','qusage'))
library(devtools)
install_github("vqv/ggbiplot")
install_github("YuLab-SMU/DOSE")
install_github("GuangchuangYu/enrichplot")
install_github("christensensm/COMPOSE")
```

