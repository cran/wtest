#' Odds ratio
#'
#' @description Function to calculate odds ratio of both main effect and pair wise interaction. Odds ratio of pair wise interaction is calculated by the coefficient of the interaction term in Logistic Regression.
#' @param data a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1).
#' @param y a numeric vector.
#' @param w.order a numeric number taking values 1 or 2. If w.order = 1, main effect is calculated. If w.order = 2, pairwise interaction effect is calculated.
#' @param which.pair a numeric vector with length = w.order indicates which pair to calculate the odds ratio.
#' @return The odds ratio of one marker
#' @export
#' @examples
#' data(mydata)
#' data(phenotype)
#' y<-as.numeric(phenotype)
#' result <- odds.ratio(mydata,y,w.order=2,which.pair=c(4,8))
#' @importFrom stats glm binomial

odds.ratio<-function(data,y,w.order,which.pair){
  if(!is.data.frame(data))
    data<-as.data.frame(data)
  if(any(is.na(data)))
    stop("NA occurs in data")
  if(!all(as.matrix(data) %in% c(0,1,2)))
    stop("all the genotypes in 'data' must be 0, 1 or 2")
  if(!is.null(which.pair) & length(which.pair)!=w.order)
    stop(gettextf("the length of 'which.pair' is %d, should equal to %d (the number of 'w.order' defined)",length(which.pair),w.order))
  if(length(y)!=nrow(data))
    stop("'data' and 'y' must have the same length")
  n.snp<-ncol(data)
  if(w.order==1){
    or.table<-table(as.matrix(data[,which.pair]),y)
    if(0 %in% or.table)
      or.table<-or.table+0.5
    if(nrow(or.table)==3)
      result.oddsratio<-(2*or.table[3,2]+or.table[2,2])*(2*or.table[1,1]+or.table[2,1])/(2*or.table[1,2]+or.table[2,2])/(2*or.table[3,1]+or.table[2,1])
    else
      result.oddsratio<-or.table[2,2]*(2*or.table[1,1]+or.table[2,1])/(2*or.table[1,2]+or.table[2,2])/or.table[2,1]
  }
  else if(w.order==2){
    result.glm<-glm(y~data[,which.pair[1]]*data[,which.pair[2]],family=binomial)$coefficients[4]
    result.oddsratio<-unname(exp(result.glm))
  }
  return(result.oddsratio)
}

