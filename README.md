EPIC package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Description
-----------

Package implementing EPIC method to estimate the proportion of immune, stromal, endothelial and cancer or other cells from bulk gene expression data. It is based on reference gene expression profiles for the main non-malignant cell types and it predicts the proportion of these cells and of the remaining "other cells" (that are mostly cancer cells) for which no reference profile is given.

This method is described in the publication from *Racle et al., 2017* available at <https://elifesciences.org/articles/26476>.

Usage
-----

The main function in this package is `EPIC`. It needs as input a matrix of the TPM (or RPKM) gene expression from the samples for which to estimate cell proportions. One can also define the reference cells to use

``` r
out <- EPIC(bulk = bulkSamplesMatrix)
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList)
```

`out` is a list containing the various mRNA and cell fractions in each samples as well as some *data.frame* of the goodness of fit.

Values of mRNA per cell and signature genes to use can also be changed:

``` r
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell = mRNA_cell_vector, sigGenes = sigGenes_vector)
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell_sub = mRNA_cell_sub_vector)
```

Installation
------------

``` r
install.packages("devtools")
devtools::install_github("GfellerLab/EPIC")
```

License
-------

EPIC can be used freely by academic groups for non-commercial purposes. The product is provided free of charge, and, therefore, on an "*as is*" basis, without warranty of any kind. Please read the file "*LICENSE*" for details.

If you plan to use EPIC (version 1.1) in any for-profit application, you are required to obtain a separate license. To do so, please contact Ece Auffarth (<eauffarth@licr.org>) at the Ludwig Institute for Cancer Research Ltd.

Contact information
-------------------

Julien Racle (<julien.racle@unil.ch>), and David Gfeller (<david.gfeller@unil.ch>).
