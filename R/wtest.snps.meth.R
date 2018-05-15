#' W-test for gene-methylation interaction analysis
#' @description {This function performs the \code{W-test} to calculate gene-methylation interactions at a (SNP,CpG) pair level for categorical
#' data sets and suitable for a genome-wide testing. This function can automatically screen and exhaustively evaluate interaction effects of
#' the SNPs and CpG sites located within a user-defined genome distance. The output can be filtered by p-values, such that only sets with smaller
#' p-value than a threshold (\code{output.pval}) will be returned.}
#' @param geno a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1). SNP names should be stored at column names of the data.
#' @param meth a data frame or matrix contains methylation data in the columns. Methylation data should be recoded as (0, 1, 2) or (0, 1). Names of CpG sites should be stored at column names of the data.
#' @param y a numeric vector composed of 0 or 1; or a factor variable with two levels.
#' @param geno.pos a data frame contains SNP names and positions in two columns.
#' @param meth.pos a data frame contains CpG names and positions in two columns.
#' @param window.size a numeric number specifies the size of genome distance. Interaction effects of the SNPs and CpG sites located within the size of genome distance will be evaluated exhaustively.
#' @param hf a data frame or matrix, contains the \emph{h} and \emph{f} values for pairwise interaction effect calculation when \emph{k} = 2 to 9. Default \emph{hf} is \emph{h} = \emph{k}/(\emph{k}-1) and \emph{f} = \emph{k}-1, where \emph{k} = 2 to 9, the first row is the \emph{h} and \emph{f} for \emph{k} = 2, and the last row is the \emph{h} and \emph{f} for \emph{k} = 9.
#' @param output.pval a p-value threshold for filtering the output. If NULL, all the results will be listed; otherwise, the function will only output the results with p-values smaller than the output.pval.
#' @param sort a logical value indicating whether or not to sort the output by p-values in ascending order. Default = TRUE.
#' @return An object \code{"wtest.snps.meth"} containing:
#'
#' \item{results}{The test results include: SNP name, CpG name, SNP position, CpG position, W value, \emph{k}, and p-value.}
#'
#' \item{hf}{The \emph{h} and \emph{f} values used for each \emph{k} in pairwise calculation, where \emph{k} = 2 to 9.}
#'
#' @details {W-test is a model-free statistical test to measure main effect or pairwise interactions in case-control studies with categorical variables.
#' Theoretically, the test statistic follows a Chi-squared distribution with \emph{f} degrees of freedom. The data-adaptive degree of freedom \emph{f},
#' and a scalar \emph{h} in the test statistics allow the W-test to correct for distributional bias due to sparse data and small sample size.
#' Let \emph{k} be the number of columns of the 2 by \emph{k} contingency table formed by a single variable or a variable pair.
#' When the sample size is large and there is no population stratification, the \emph{h} and \emph{f} will approximate well to the theoretical
#' value \emph{h} = (\emph{k}-1)/\emph{k}, and \emph{f} = \emph{k}-1. When sample size is small and there is population stratification, the \emph{h} and
#' \emph{f} will vary to correct for distributional bias caused by the data structure.}
#'
#' @examples
#' data(SNP_pos)
#' data(CpG_pos)
#' data(genotype)
#' data(methylation)
#' data(phenotype2)
#'
#' w <- 13000
#'
#' # Recode methylation data
#' methylation <- methylation.recode(methylation)
#'
#' ## Step 1. HF Calculation.
#' # Please note that parameter B is recommended to be greater than 400.
#' hf.pair <- hf.snps.meth(B = 80, geno = genotype, meth = methylation, y = phenotype2,
#'                         geno.pos = SNP_pos, meth.pos = CpG_pos, window.size = w)
#'
#' ## Step 2. Application
#' result <- wtest.snps.meth(geno = genotype, meth = methylation, y = phenotype2, geno.pos = SNP_pos,
#'                           meth.pos = CpG_pos, window.size = w, hf = hf.pair, output.pval = 0.1)
#'
#' @export
#' @author Rui Sun, Maggie Haitian Wang
#' @references Maggie Haitian Wang, Rui Sun, Junfeng Guo, Haoyi Weng, Jack Lee, Inchi Hu, Pak Sham and Benny C.Y. Zee (2016). A fast and powerful W-test for pairwise epistasis testing. Nucleic Acids Research.doi:10.1093/nar/gkw347.
#' @seealso \code{\link{wtest}}, \code{\link{hf.snps.meth}}
#' @importFrom utils combn
#' @importFrom stats pchisq

wtest.snps.meth <- function(geno, meth, y, geno.pos, meth.pos, window.size = 1e4, hf = "default.hf",
                    output.pval = NULL, sort = TRUE){
  suppressWarnings(if(typeof(hf) == "character"){hf = array(c(0.5,0.667,0.75,0.8,0.833,0.857,0.875,0.889,1:8), dim = c(8,2))}else{hf = hf[,2:3]})
  if(is.data.frame(geno))
    geno <- as.matrix(geno)
  if(any(is.na(geno)))
    stop("NA occurs in data.genotype")
  if(is.data.frame(meth))
    meth <- as.matrix(meth)
  if(any(is.na(meth)))
    stop("NA occurs in data.methylation")
  if(!all(geno %in% c(0,1,2)))
    stop("all the genotypes in 'data.genotype' must be 0, 1 or 2")
  if(length(y)!=nrow(geno) || nrow(meth)!=nrow(geno))
    stop("'data.genotype', 'data.methylation' and 'y' must have the same length")

  snp.names <- colnames(geno)
  cpg.names <- colnames(meth)
  l1 <- match(snp.names, geno.pos[,1])
  l2 <- match(cpg.names, meth.pos[,1])
  if(any(is.na(l1)))
    stop("missing SNP position exists")
  if(any(is.na(l2)))
    stop("missing CpG position exists")
  geno.pos <- geno.pos[l1,]
  meth.pos <- meth.pos[l2,]

  index.set<-data.frame()
  for(i in 1:nrow(geno.pos)){
    index <- which(abs(geno.pos[,2][i] - meth.pos[,2]) <= window.size)
    if(length(index)){
      index.i <- cbind(i, index)
      index.set <- rbind(index.set, index.i)
    }
  }
  set <- apply(index.set, 1, list)

  n.snp <- ncol(geno)
  n.cpg <- ncol(meth)
  result <- lapply(set, x2.set, geno, meth, y)
  result.all <- do.call(rbind,result)
  x2.column <- 3
  df.column <- x2.column+1
  pval.column <- x2.column+2
  w.value <- result.all[,x2.column] * hf[result.all[,df.column],1]
  p.value <- pchisq(w.value, df = hf[result.all[,df.column],2], lower.tail = F)
  result.all[,x2.column] <- w.value
  result.all <- cbind(result.all,p.value)
  k <- result.all[,df.column]+1
  result.all[,df.column] <- k
  result.all <- as.data.frame(result.all)

  colnames(result.all) <- c("SNP","CpG","w","k","p-value")
  result.all[,1] <- snp.names[result.all[,1]]
  result.all[,2] <- cpg.names[result.all[,2]]
  if(!is.null(output.pval)){
    l.output.pval <- which(result.all[,pval.column] < output.pval)
    result.all <- result.all[l.output.pval,]
  }
  if(sort){
    l.order <- order(result.all[,pval.column], decreasing=F)
    result.all <- result.all[l.order,]
  }
  result.all$SNP.BP<-geno.pos[match(result.all[,1],geno.pos[,1]),2]
  result.all$CpG.BP<-meth.pos[match(result.all[,2],meth.pos[,1]),2]
  result.all<-result.all[,c(1,2,6,7,3,4,5)]
  k <- c(2:(nrow(hf)+1))
  hf <- cbind(k,hf)
  colnames(hf)<-c("k","h","f")
  return(list(results = result.all, hf = hf))
}
