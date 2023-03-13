## ---- eval = FALSE------------------------------------------------------------
#  # library(EPIC) ## If the package isn't loaded (or use EPIC::EPIC and so on).
#  out <- EPIC(bulk = bulkSamplesMatrix)
#  out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList)

## ---- eval = FALSE------------------------------------------------------------
#  out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell = mRNA_cell_vector, sigGenes = sigGenes_vector)
#  out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell_sub = mRNA_cell_sub_vector)

## ---- eval = FALSE------------------------------------------------------------
#  ?EPIC::EPIC
#  ?EPIC::EPIC.package

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

