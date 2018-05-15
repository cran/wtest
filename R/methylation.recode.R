#' Recode methylation data
#' @description {This function performs redocing for methylation data by a K-means clustering method. The methylation data set will be recoded into two (high and low) levels according to the methylated level of each CpG site.}
#' @param data a data frame or matrix contains methylation data in the columns.
#' @examples
#' data(methylation)
#' data.recoded <- methylation.recode(methylation)
#' @export
#' @importFrom stats kmeans
#'

methylation.recode<-function(data){
  data.recode<-apply(data,2,cluster)
  return(data.recode)
}

