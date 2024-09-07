#' TCGA BRCA Expression Data
#'
#' A dataset containing the top 2000 features of mRNA, miRNA, and DNAm from the TCGA BRCA data.
#'
#' @format A list with three data frames for mRNA, miRNA, and DNAm features:
#' \describe{
#'   \item{gene}{A data frame with mRNA features}
#'   \item{miRNA}{A data frame with miRNA features}
#'   \item{methy}{A data frame with DNAm features}
#' }
#' @source \url{https://www.cancer.gov/tcga}
"tcga_brca"

#' TCGA BRCA Clinical Data
#'
#' A dataset containing clinical data from the TCGA BRCA study.
#'
#' @format A data frame with clinical data.
#' @source \url{https://www.cancer.gov/tcga}
"tcga_brca_clinical"
