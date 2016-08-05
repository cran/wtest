#' W-test
#'
#' @description {This function performs the \code{W-test} to calculate main effect or pairwise interactions in case-control studies
#' for categorical data sets. The test measures target variables' distributional difference between cases and controls via a combined
#' log of odds ratio. It follows a chi-squared probability distribution with data-adaptive degrees of freedom. For pairwise interaction
#' calculation, the user has 3 options: (1) calculate a single pair's W-value, (2) calculate pairwise interaction for a list of variables,
#' which p-values are smaller than a threshold (\code{input.pval}); (3) calculate the pairwise interaction exhaustively for all variables.
#' For both main effect and interaction effect calculation, the output can be filtered by p-values, such that only sets with smaller p-value
#' than a threshold (\code{output.pval}) will be returned.}
#' @param data a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1).
#' @param y a numeric vector composed of 0 or 1; or a factor variable with two levels.
#' @param w.order an integer value of 0 or 1. \code{w.order} = 1 for main effect calculation; \code{w.order} = 2, for pairwise calculation.
#' @param hf1 a data frame or matrix, contains the h and f values for main effect (\code{w.order} =1) calculation at the number of categorical combinations (k) = 2 or 3. \code{Default.hf1} = array(c(0.5, 0.667, 1, 2), dim=c(2,2)), in which the first row is the h and f for k=2, and second row is the h and f for k=3.
#' @param hf2 a data frame or matrix, contains the h and f values for pairwise interaction effect calculation (\code{w.order}=2) when k = 2 to 9. \code{Default.hf2} = array(c(0.5, 0.667, 0.75, 0.8, 0.833, 0.857, 0.875, 0.889, 1:8), dim=c(8,2)), the first row is the h and f for k=2, and the last row is the h and f for k=9.
#' @param which.pair a numeric vector, with length = \code{w.order}. It contains the column number of the variable set to calculate. If which.pair is specified, the w.value for that set is returned. Default \code{which.pair} = NULL, when main or interaction effect will be calculated exhaustively.
#' @param output.pval a p-value threshold for filtering the output. If NULL, all the results will be listed; otherwise, the function will only output the results with p-values smaller than the output.pval.
#' @param sort a logical value indicating whether or not to sort the output by p-values in ascending order. Default = TRUE.
#' @param input.pval a p-value threshold to select markers for pairwise calculation, used only when \code{w.order} = 2. When specified, only markers with main effect p-value smaller than input.pval will be passed to interaction effect calculation. Default = 0.10. Set \code{input.pval} = NULL or 1 for exhaustive pairwise calculation.
#' @param input.poolsize an integer, with value less than the number of input variables. It is an optional filter to control the maximum number of variables to execute pairwise calculation, used only when \code{w.order} = 2. It selects top input.poolsize number of variables to calculate pairwise interactions. It can be used separately or jointly with \code{input.pval}, whichever gives smaller input variable pool size. Default = 50. Set \code{input.poolsize} = NULL for exhaustive pairwise calculation. It can be useful when the user is exploring the data, and there may be a large number of variables with extremely small main effect p-values.
#' @return An object \code{"wtest"} containing:
#'
#' \item{order}{the "w.order" specified.}
#'
#' \item{results}{When order = 1, the test results include: the ID of SNP, the W value, k, and p value. When order = 2 and which.pair = NULL, the test results include: (for pair, column 1-3) [ pair name, W-value, k, p-value]; (for first variable in the pair, column 4-6) [W-value, k, p-value ]; (for second variable in the pair, column 7-8) [W-value, k, p-value].}
#'
#' \item{hf1}{The h and f values used in main effect calculation.}
#'
#' \item{hf2}{The h and f values used in pairwise calculation.}
#'
#' @details {W-test is a model-free statistical test to measure main effect or pairwise interactions in case-control studies with categorical variables.
#' Theoretically, the test statistic follows a Chi-squared distribution with f degrees of freedom. The data-adaptive degree of freedom f,
#' and a scalar h in the test statistics allow the W-test to correct for distributional bias due to sparse data and small sample size.
#' Let k be the number of columns of the 2 by k contingency table formed by a single variable or a variable pair.
#' When the sample size is large and there is no population stratification, the h and f will approximate well to the theoretical
#' value h = (k-1)/k, and f = k-1. When sample size is small and there is population stratification, the h and
#' f will vary to correct for distributional bias caused by the data structure.}
#'
#' {When \code{w.order} =2, the \code{wtest()} will automatically calculate the main effect first and then do a pre-filter before calculating interactions.
#' This filtering is to avoid overloading the memory before having a better understanding of the data. User can specify a smaller input.pval such as 0.05 or 0.001
#' for less output, or \code{input.pval}=1 or NULL for exhaustive pairwise calculation. Another optional filter is input.poolsize. It will take the top \code{input.poolsize}
#' number of variables to calculated pairwise effect exhaustively, selected by smallest p-value; when used together with \code{input.pval}, the smaller set will be passed to pairwise calculation.}
#'
#' @examples
#' data(mydata)
#' data(phenotype)
#' hf1<-hf.calculation(data = mydata, w.order = 1, B = 100)
#' hf2<-hf.calculation(data = mydata, w.order = 2, B = 50)
#' w1<-wtest(mydata, phenotype, w.order=1, hf1 = hf1)
#' w2<-wtest(mydata, phenotype, w.order=2, input.pval = 0.5, output.pval = 0.01, hf1 = hf1, hf2 = hf2)
#' w.pair<-wtest(mydata, phenotype, w.order=2, which.pair=c(10,13), hf2 = hf2)
#' @export
#' @author Rui Sun, Maggie Haitian Wang
#' @references Maggie Haitian Wang, Rui Sun, Junfeng Guo, Haoyi Weng, Jack Lee, Inchi Hu, Pak Sham and Benny C.Y. Zee (2016). A fast and powerful W-test for pairwise epistasis testing. Nucleic Acids Research.doi:10.1093/nar/gkw347.
#' @seealso \code{\link{hf.calculation}}, \code{\link{w.diagnosis}}, \code{\link{w.qqplot}}
#' @importFrom utils combn
#' @importFrom stats pchisq

wtest<-function(data,y,w.order=c(1,2),hf1="default.hf1",hf2="default.hf2",
                which.pair=NULL,output.pval=NULL,sort=TRUE,input.pval=0.10,input.poolsize=50){
  suppressWarnings(if(hf1=="default.hf1"){hf1=array(c(0.5,0.667,1,2),dim=c(2,2))}else{hf1=hf1[,2:3]})
  suppressWarnings(if(hf2=="default.hf2"){hf2=array(c(0.5,0.667,0.75,0.8,0.833,0.857,0.875,0.889,1:8),dim=c(8,2))}else{hf2=hf2[,2:3]})
  if(is.data.frame(data))
    data<-as.matrix(data)
  if(any(is.na(data)))
    stop("NA occurs in data")
  if(!all(data %in% c(0,1,2)))
    stop("all the genotypes in 'data' must be 0, 1 or 2")
  if(!is.null(which.pair) & length(which.pair)!=w.order)
    stop(gettextf("the length of 'which.pair' is %d, should equal to %d (the number of 'w.order' defined)",length(which.pair),w.order))
  if(length(y)!=nrow(data))
    stop("'data' and 'y' must have the same length")
  n.snp<-ncol(data)
  if(!is.null(which.pair)){
    set<-list(which.pair)
  }else if(w.order==1){
    set<-lapply(1:n.snp,list)
  }else if(w.order==2){
    if(is.null(input.pval) & is.null(input.poolsize)){
      set<-apply(t(combn(n.snp,2)),1,list)
    }else{
      input.pval<-ifelse(is.null(input.pval),1,input.pval)
      input.poolsize<-ifelse(is.null(input.poolsize),n.snp,input.poolsize)
      set.order1<-lapply(1:n.snp,list)
      result.order1<-lapply(set.order1,x2,data,y,1)
      result.order1.all<-do.call(rbind,result.order1)
      w.value.order1<-result.order1.all[,2]*hf1[result.order1.all[,3],1]
      p.value.order1<-pchisq(w.value.order1,df=hf1[result.order1.all[,3],2],lower.tail=F)
      result.order1.all[,2]<-w.value.order1
      result.order1.all[,3]<-result.order1.all[,3]+1
      result.order1.all<-cbind(result.order1.all,p.value.order1)
      l.select<-which(result.order1.all[,4]<input.pval)
      if(length(l.select)>input.poolsize){
        result.order1.rank<-result.order1.all[order(result.order1.all[,4],decreasing=F),]
        l.select<-result.order1.rank[1:input.poolsize,1]
      }
      set<-apply(t(combn(l.select,2)),1,list)
    }
  }
  result<-lapply(set,x2,data,y,w.order)
  result.all<-do.call(rbind,result)
  x2.column<-ifelse(w.order==1,2,3)
  if(w.order==1){
    hf<-hf1
  }else {
    hf<-hf2
  }
  df.column<-x2.column+1
  pval.column<-x2.column+2
  w.value<-result.all[,x2.column]*hf[result.all[,df.column],1]
  p.value<-pchisq(w.value,df=hf[result.all[,df.column],2],lower.tail=F)
  result.all[,x2.column]<-w.value
  result.all<-cbind(result.all,p.value)
  k<-result.all[,df.column]+1
  result.all[,df.column]<-k
  marker.names<-colnames(data)
  result.all<-as.data.frame(result.all)
  if(w.order==2){
    colnames(result.all)<-c("marker1","marker2","w","k","p-value")
    if(!(is.null(input.pval) & is.null(input.poolsize)) & is.null(which.pair)){
      result.all<-cbind(result.all,result.order1.all[result.all[,1],2:4],result.order1.all[result.all[,2],2:4])
      colnames(result.all)<-c("marker1","marker2","w","k","pair.p-value","marker1.w","marker1.k","marker1.p-value","marker2.w","marker2.k","marker2.p-value")
    }
    result.all[,1]<-marker.names[result.all[,1]]
    result.all[,2]<-marker.names[result.all[,2]]
  }else{
    result.all[,1]<-marker.names[result.all[,1]]
    colnames(result.all)<-c("marker","w","k","p-value")
  }
  if(!is.null(output.pval)){
    l.output.pval<-which(result.all[,pval.column]<output.pval)
    result.all<-result.all[l.output.pval,]
  }
  if(sort){
    l.order<-order(result.all[,pval.column],decreasing=F)
    result.all<-result.all[l.order,]
  }
  k<-c(2:(nrow(hf1)+1))
  hf1<-cbind(k,hf1)
  k<-c(2:(nrow(hf2)+1))
  hf2<-cbind(k,hf2)
  return(list(order=w.order, results = result.all,hf1=hf1,hf2=hf2))
}
