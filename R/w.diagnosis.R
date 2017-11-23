w.null<-function(index,w.order,n.sample,n.marker,data,hf){
  if(n.sample < nrow(data)){
    sample.select<-sample(nrow(data),n.sample)
    data<-data[sample.select,]
  }
  if(n.marker < ncol(data)){
    snp.select<-sample(ncol(data),n.marker)
    data<-data[,snp.select]
  }
  y<-sample(0:1,n.sample,replace=T)
  w.order<-unlist(w.order)
  set<-apply(t(combn(n.marker,w.order)),1,list)
  result<-lapply(set,x2,data,y,w.order)
  result.all<-do.call(rbind,result)
  df.column<-ifelse(w.order==1,3,4)
  w.column<-df.column-1
  w.value<-result.all[,w.column]*hf[result.all[,df.column],1]
  df.k<-hf[result.all[,df.column],2]
  w.null<-cbind(w.value,df.k)
  return(w.null)
}

#' Plots for w values diagnosis
#'
#' @description Function to plot the distribution of w values against its theoretical distribution for W-test diagnosis.
#' @param data a data frame or matrix contains genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1).
#' @param w.order an integer value of 0 or 1. w.order = 1 for main effect calculation; w.order = 2, for pairwise calculation.
#' @param n.rep a numeric value, the number of replicates.
#' @param n.sample a numeric value, the number of samples to use in bootstrapping. Default is the total number of samples in the data.
#' @param n.marker a numeric value, the number of markers to use in bootstrapping. Default is the total number of markers.
#' @param hf1 a data frame or matrix, contains the h and f values for main effect (w.order =1) calculation at the number of categorical combinations (k) = 2 or 3. Default.hf1= array(c(0.5, 0.667, 1, 2), dim=c(2,2)), in which the first row is the h and f for k=2, and second row is the h and f for k=3.
#' @param hf2 a data frame or matrix, contains the h and f values for pairwise interaction effect calculation (w.order=2) when k = 2 to 9. Default.hf2 = array(c(0.5, 0.667, 0.75, 0.8, 0.833, 0.857, 0.875, 0.889, 1:8), dim=c(8,2)), the first row is the h and f for k=2, and the last row is the h and f for k=9.
#' @param ... graphical parameters.
#'
#' @details {This function evaluates the input w values of main or interaction effects using a set of null Y by the \code{W-test}, and plots the w values probabilities against its theoretical Chi-squared distribution density with f degrees of freedom indexed by k.
#' Since h and f are dataset adaptive, estimated hf is strongly recommended.}
#'
#' @examples
#' data(mydata)
#' # Please note that parameter B is recommended to be greater than 400.
#' hf1<-hf.calculation(data = mydata, w.order = 1, B = 100)
#' hf2<-hf.calculation(data = mydata, w.order = 2, B = 50)
#' w.diagnosis(mydata, w.order = 1, n.rep = 100, hf1 = hf1, main=NULL, xlab=NULL, ylab=NULL)
#' w.diagnosis(mydata, w.order = 2, n.rep = 100, hf2 = hf2, main=NULL, xlab=NULL, ylab=NULL)
#' @export
#' @author Rui Sun, Maggie Haitian Wang
#' @references Maggie Haitian Wang, Rui Sun, Junfeng Guo, Haoyi Weng, Jack Lee, Inchi Hu, Pak Sham and Benny C.Y. Zee (2016). A fast and powerful W-test for pairwise epistasis testing. Nucleic Acids Research.doi:10.1093/nar/gkw347.
#' @seealso \code{\link{wtest}}, \code{\link{hf.calculation}}, \code{\link{w.qqplot}}
#' @importFrom graphics par hist lines text mtext plot legend
#' @importFrom stats rchisq density
#' @importFrom utils combn

w.diagnosis<-function(data, w.order=c(1,2), n.rep=10, n.sample=nrow(data), n.marker=ncol(data),
                 hf1="default.hf1", hf2="default.hf2", ...){
  suppressWarnings(if(hf1=="default.hf1"){hf1=array(c(0.5,0.667,1,2),dim=c(2,2))}else{hf1=hf1[,2:3]})
  suppressWarnings(if(hf2=="default.hf2"){hf2=array(c(0.5,0.667,0.75,0.8,0.833,0.857,0.875,0.889,1:8),dim=c(8,2))}else{hf2=hf2[,2:3]})
  if(is.data.frame(data))
    data<-as.matrix(data)
  if(any(is.na(data)))
    stop("NA occurs in data")
  if(!all(data %in% c(0,1,2)))
    stop("all the genotypes in 'data' must be 0, 1 or 2")
  if(w.order==1){
    hf<-hf1
  }else {
    hf<-hf2
  }
  index<-apply(array(1:n.rep,dim=c(n.rep,1)),1,list)
  w.value<-lapply(index,w.null,w.order,n.sample,n.marker,data,hf)
  w.all<-do.call(rbind,w.value)
  k.all<-match(unique(w.all[,2]),hf[,2])+1
  n.p<-length(k.all)
  op<-par(no.readonly=TRUE)
  par(mfrow=c(ceiling(n.p/2),min(2,n.p)),oma = c(3, 2, 2, 1),mai=c(0.4,0.4,0.1,0.1),xpd=T)
  for(k in sort(k.all,decreasing=T)){
    o.v<-w.all[which(w.all[,2]==hf[k-1,2]),1]
    e.v<-rchisq(length(o.v),df=hf[k-1,2])
    hist(o.v,freq=F,ylim=c(0,max(max(density(o.v,na.rm=T)[[2]]),max(density(e.v)[[2]]))),...)
    lines(density(e.v,na.rm=T),lty=1,col="red")
    text(x=max(o.v)*4/5,y=max(density(o.v,na.rm=T)[[2]])*4/5,paste0("k = ",k,"; df = ",round(hf[k-1,2],2)))
  }
  mtext("w value",outer=T,side=1,cex=0.7)
  mtext("density",outer=T,side=2,cex=0.7)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend=c("expected ","observed"), ncol=2, cex=1,
         inset = c(0, 0), bty="n", pch=c(NA,0), col=c("red","black"),lwd=c(1.5,NA),lty=c(1,NA))
  par(op)
}
