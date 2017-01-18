#' EPIC: a package to Estimate the Proportion of Immune and Cancer cells from
#'  tumor gene expression data.
#'
#' EPIC package provides the function and immune cell reference profiles to
#' estimate the proportion of immune and other cells from bulk gene expression
#' data.
#'
#' @section EPIC functions:
#' \code{\link{EPIC}} is the main function to call to estimate the
#'  various cells proportions from a bulk sample.
#'
#' @section Included datasets:
#' \code{\link{BRef}}, \code{\link{BRef.tpm}}: reference profiles from
#'    circulating immune cells.
#'
#'  \code{\link{TRef.tpm}}: reference profiles from immune cells obtained
#'  from single cell data of tumor infiltrating cells from melanoma patients.
#'
#' \code{\link{hoek_data}}: example dataset containing data from Hoek et al,
#'  2015, PLoS One.
#'
#' \code{\link{mRNA_cell_default}}: values of mRNA per cell for the main cell
#'  types.
#'
#' @section Authors:
#' Julien Racle <\email{julien.racle@unil.ch}> and David Gfeller <\email{
#'  david.gfeller@unil.ch}>.
#'
#' @docType package
#' @name EPIC.package
NULL


#' Reference profiles from circulating immune cells.
#'
#' A dataset containing the reference profiles obtained from immune cell
#' samples of \emph{B-}, \emph{NK-}, \emph{T-cells}, \emph{Monocytes}
#' and \emph{Neutrophils} purified from PBMC or whole blood.
#'
#' The original samples were obtained from healthy donors and donors after
#' influenza vaccination or with diabetes, sepsis or multiple sclerosis.
#'
#' @section Similar datasets:
#'  \code{BRef} (main reference profiles, using data from sources 1-3 below)
#'
#'
#' @format A list of 3 elements: \describe{ \item{$refProfiles,
#'   $refProfiles.var}{Matrices (nGenes x nRefCells) of the normalized gene
#'   expression from the reference cells and the variability of this gene
#'   expression for each gene and each cell type} \item{$sigGenes}{A list of 100
#'   signature genes used to deconvolve the cell proportions} }
#'
#' @source \enumerate{
#'  \item \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-64655},
#'  \item \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-60424/},
#'  \item \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51984}
#'  }
"BRef"

#' @section Similar datasets:
#' \code{BRef.tpm} (reference profiles based on same data as \code{BRef},
#'  but given in TPM counts instead of raw counts)
#' @rdname BRef
"BRef.tpm"

#' Reference profiles obtained from single cell data of tumor infiltrating
#' cells.
#'
#' A dataset containing the reference profiles given in TPM from various cell
#' types: \emph{B-}, \emph{NK-}, \emph{T-cells} and \emph{Macrophages}.
#'
#' These were obtained from single-cell RNA-seq data from 9 donors from
#' the publication of Tirosh et al., 2016, Science. The samples
#' come from cancer metastases of melanoma (extracted from primary tumors
#' and non-lymphoid tissue metastases). The classification for each sample with
#' respect to each cell type is the one given by Tirosh et al.
#'
#' @format A list of 3 elements: \describe{ \item{$refProfiles,
#'   $refProfiles.var}{Matrices (nGenes x nRefCells) of the normalized gene
#'   expression from the reference cells and the variability of this gene
#'   expression for each gene and each cell type} \item{$sigGenes}{A list of 80
#'   signature genes used to deconvolve the cell proportions} }
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}
"TRef.tpm"

#' Values of mRNA / cell for the main cell types.
#'
#' These values have been obtained from experiments (see \cite{EPIC} publication).
#' For the other uncharacterized cells, we use a value of 0.4 as described
#' in \cite{EPIC} publication. For macrophages we don't have specific values but
#' assumed here it is the same value as for monocytes.
#'
#' @format A named numeric vector of the relative amount of mRNA per cell type.
#'  There are two additional "special cell types": the \emph{otherCells} which
#'  correspond to the uncharacterized cells present in the sample but without
#'  any reference profile and the \emph{default} which is the default value used
#'  for cells with reference profiles but without a value specified in the
#'  \code{mRNA_cell_default} vector.
"mRNA_cell_default"

#' Example dataset containing data from Hoek et al, 2015, PLoS One.
#'
#' This dataset contains a subset of the full Hoek et al data. It contains only
#' the data from the two healthy donors PBMC before influenza vaccination.
#'
#' @format This is a list of 3 elements: \describe{
#'  \item{$rawCounts}{(matrix of 51574 genes x 2 donors) The raw read counts
#'    from the two donors. It has been obtained by mapping the original
#'    fastq files to \emph{hg19} genome with help of \emph{tophat} and
#'    \emph{htseq-count}.}
#'  \item{$cellFractions.obs}{(matrix of 2 donors x 6 cell types) The
#'    proportions of the different immune cells in the PBMC from the 2 donors,
#'    as measured by FACS by Hoek et al.}
#'  \item{$cellFractions.pred}{(matrix of 2 donors x 7 cell types) The
#'    proportions of the different immune cells and of a potential additional
#'    uncharacterized cell, as predicted by EPIC based on the rawCounts and
#'    the profiles \code{reference=BRef}.}
#' }
#'
#' @source The description of this data can be found here:
#'  \href{http://dx.doi.org/10.1371/journal.pone.0118528}{link to paper}
#'  and \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-64655}{link
#'  to data}.
"hoek_data"

