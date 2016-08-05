#' plots for w p-values diagnosis
#' @description To draw a Q-Q plot for W-test
#' @param data a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1).
#' @param w.order a numeric number taking values 1 or 2. If w.order = 1, main effect is calculated. If w.order = 2, pairwise interaction effect is calculated.
#' @param hf1 a data frame or matrix, specifying the h and f parameters for order 1 calculation. Default values = HF1= array(c(0.5,0.667,1,2), dim=c(2,2)).
#' @param hf2 a data frame or matrix, specifying the h and f parameters for order 2 calculation. Default values = HF2 = array(c(0.5,0.667,0.75,0.8,0.833,0.857,0.875,0.889,1:8), dim=c(8,2)).
#' @param input.poolsize a numeric number; The maximum number of values used to calculate first or second order effects. The default is 200.
#' @param ... graphical parameters.
#' @return Q-Q plot
#' @details
#' The Q-Q plot for W-test is to use a set of randomly generated y as phenotype, to test its null distribution compared with chi-square distribution. To fit different type of data and adjust the distribution, h and f parameters are strongly recommended to calculate instead of the default hf1 and hf2 for first and second order, respectively.
#'
#' The input.poolsize is suggested to set as 1000 for w.order = 1 and 200 for w.order = 2.
#'
#' @examples
#' data(mydata)
#' w.qqplot(mydata,w.order=2,cex=.5)
#' abline(0,1)
#' @export
#' @importFrom utils combn
#' @importFrom stats pchisq runif qqplot
#'
w.qqplot<-function(data, w.order=c(1,2), input.poolsize=50, hf1="default.hf1", hf2="default.hf2",
                   ...){
  suppressWarnings(if(hf1=="default.hf1"){hf1=array(c(0.5,0.667,1,2),dim=c(2,2))}else{hf1=hf1[,2:3]})
  suppressWarnings(if(hf2=="default.hf2"){hf2=array(c(0.5,0.667,0.75,0.8,0.833,0.857,0.875,0.889,1:8),dim=c(8,2))}else{hf2=hf2[,2:3]})
  if(is.data.frame(data))
    data<-as.matrix(data)
  if(any(is.na(data)))
    stop("NA occurs in data")
  if(!all(data %in% c(0,1,2)))
    stop("all the genotypes in 'data' must be 0, 1 or 2")
  y<-sample(0:1,nrow(data),replace=T)
  n.snp<-ncol(data)
  if(w.order==1){
    set<-lapply(1:n.snp,list)
  }else if(w.order==2){
    if(n.snp<=input.poolsize){
      set<-apply(t(combn(n.snp,2)),1,list)
    }else{
      l.select<-sample(1:n.snp,input.poolsize,replace=F)
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
  p.value.observed<-pchisq(w.value,df=hf[result.all[,df.column],2],lower.tail=F)
  p.value.expected<-runif(length(p.value.observed))
  qqplot(-log10(p.value.expected),-log10(p.value.observed),...)
}
