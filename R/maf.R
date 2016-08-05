#' Minor allele frequency
#'
#' @description Function to calculate minor allele frequencies.
#' @param data a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1).
#' @param which.marker a numeric vector with length = w.order indicates which SNP to calculate the MAF. When which.marker = NULL, MAF of all the markers is calculated. Default is NULL.
#' @return The MAF of one marker
#' @examples
#' data(mydata)
#' result<-maf(mydata, which.marker=10)
#' @export
maf<-function(data,which.marker=NULL){
  if(!is.null(which.marker) & !is.numeric(which.marker))
    stop("the 'which.marker' should be numeric!")
  if(!is.matrix(data))
    data<-as.matrix(data)
  if(any(is.na(data)))
    stop("NA occurs in data")
  if(!all(data %in% c(0,1,2)))
    stop("all the genotypes in data should be 0, 1 or 2")
  if(!is.null(which.marker)){
    result.maf<-mean(data[,which.marker])/2
  }else{
    result.maf<-colMeans(data)/2
  }
  return(result.maf)
}
