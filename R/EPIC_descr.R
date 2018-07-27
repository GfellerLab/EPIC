#' EPIC: a package to Estimate the Proportion of Immune and Cancer cells from
#'  tumor gene expression data.
#'
#' EPIC package provides the function and cell reference profiles to
#' estimate the proportion of immune, stromal, endothelial and cancer or other
#' cells from bulk gene expression data.
#'
#' See the package \link[=../doc/info.html]{vignette} and function definitions
#' below.
#'
#' @section EPIC functions:
#' \code{\link{EPIC}} is the main function to call to estimate the
#'  various cells proportions from a bulk sample.
#'
#' @section Included datasets:
#' \code{\link{BRef}}: reference profiles from circulating immune cells.
#'
#'  \code{\link{TRef}}: reference profiles from tumor infiltrating non-malignant
#'  cells obtained from single cell data of melanoma patients.
#'
#' \code{\link{melanoma_data}}: example dataset containing data from lymph nodes
#'  from patients with metastatic melanoma.
#'
#' \code{\link{mRNA_cell_default}}: values of mRNA per cell for the main cell
#'  types.
#'
#' @section References:
#' Racle, J., Jonge, K. de, Baumgaertner, P., Speiser, D.E., and
#' Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell
#' types from bulk tumor gene expression data. \emph{eLife Sciences, 6}, e26476
#' (\url{https://elifesciences.org/articles/26476})
#'
#' @section Authors:
#' Julien Racle <\email{julien.racle@unil.ch}> and David Gfeller
#' <\email{david.gfeller@unil.ch}>.
#'
#' @docType package
#' @name EPIC.package
NULL


#' Reference profiles from circulating immune cells.
#'
#' A dataset containing the reference profiles obtained from immune cell
#' samples of \emph{B cells}, \emph{CD4 T cells}, \emph{CD8 T cells},
#' \emph{Monocytes}, \emph{NK cells} and \emph{Neutrophils}, purified from
#' PBMC or whole blood.
#'
#' The original samples were obtained from healthy donors and donors after
#' influenza vaccination or with diabetes, sepsis or multiple sclerosis.
#'
#' @format A list of 3 elements: \describe{ \item{$refProfiles,
#'   $refProfiles.var}{Matrices (nGenes x nRefCells) of the gene expression (in
#'   TPM counts) from the reference cells and the variability of this gene
#'   expression for each gene and each cell type} \item{$sigGenes}{A list of
#'   signature genes used to deconvolve the cell proportions} }
#'
#' @source \enumerate{
#'  \item \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-64655},
#'  \item \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-60424/},
#'  \item \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51984}
#'  }
"BRef"


#' Reference profiles obtained from single cell data of tumor infiltrating
#' cells.
#'
#' A dataset containing the reference profiles obtained from various
#' tumor infiltrating non-malignant cell types: \emph{B cells},
#' \emph{cancer-associated fibroblasts}, \emph{CD4 T cells}, \emph{CD8 T cells},
#' \emph{endothelial cells}, \emph{Macrophages} and \emph{NK cells}.
#'
#' These were obtained from single-cell RNA-seq data from 9 donors from
#' the publication of \href{http://science.sciencemag.org/content/352/6282/189.full}{
#' \cite{Tirosh et al., 2016, Science}}. The samples
#' come from melanoma tumors (extracted from primary tumors and non-lymphoid
#' tissue metastases). The classification for each sample with
#' respect to each cell type is the one given by Tirosh et al., except for
#' the CD4 T cells and CD8 T cells, that were identified from the T cells based
#' on the expression of CD4, CD8A and CD8B as described in
#' \href{https://elifesciences.org/articles/26476}{
#'  \cite{Racle et al., 2017, eLife}}.
#'
#' @format A list of 3 elements: \describe{ \item{$refProfiles,
#'   $refProfiles.var}{Matrices (nGenes x nRefCells) of the gene expression (in
#'   TPM counts) from the reference cells and the variability of this gene
#'   expression for each gene and each cell type} \item{$sigGenes}{A list of
#'   signature genes used to deconvolve the cell proportions} }
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}
"TRef"

#' Values of mRNA / cell for the main cell types.
#'
#' These values have been obtained from experiments (see
#' \href{https://elifesciences.org/articles/26476}{
#'  \cite{Racle et al., 2017, eLife}}).
#' For the other uncharacterized cells, we use a value of 0.4 as described
#' in \href{https://elifesciences.org/articles/26476}{
#'  \cite{Racle et al., 2017, eLife}}.
#' For macrophages we don't have specific values but assumed here it is the
#' same value as for monocytes.
#'
#' @format A named numeric vector of the relative amount of mRNA per cell type.
#'  There are two additional "special cell types": the \emph{otherCells} which
#'  correspond to the uncharacterized cells present in the sample but without
#'  any reference profile and the \emph{default} which is the default value used
#'  for cells with reference profiles but without a value specified in the
#'  \code{mRNA_cell_default} vector.
#'
#' @source \url{https://elifesciences.org/articles/26476}
"mRNA_cell_default"

#' Example dataset containing data from lymph nodes from patients with
#' metastatic melanoma.
#'
#' This is the dataset measured in \href{https://elifesciences.org/articles/26476}{
#'  \cite{Racle et al., 2017, eLife}}. It contains
#' the gene expression from lymph node samples from four patients with melanoma,
#' and it contains also the proportions of the main immune cell types and of
#' melanoma cells, as measured by FACS.
#'
#' @format This is a list of 3 elements: \describe{
#'  \item{$counts}{(matrix of 49902 genes x 4 donors) The TPM normalized counts
#'    from the four donors. It has been obtained by mapping RNA-seq data to
#'    \emph{hg19} genome with help of \emph{RSEM}. Ensembl ID were then
#'    converted to gene names, and genes with duplicated entries were
#'    merged together by summing their counts.}
#'  \item{$cellFractions.obs}{(matrix of 4 donors x 6 cell types) The
#'    proportions of the different cell types measured by FACS (the
#'    "other_cells" correspond to the live cells without any marker of the
#'    other given cell types).}
#'  \item{$cellFractions.pred}{(matrix of 4 donors x 8 cell types) The
#'    proportions of the different cell types, as predicted by EPIC based on
#'    the reference profiles \code{TRef}.}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93722}
"melanoma_data"

