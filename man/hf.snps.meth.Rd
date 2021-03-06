% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hf.snps.meth.R
\name{hf.snps.meth}
\alias{hf.snps.meth}
\title{Parameter Estimation for W-test Probability Distribution in Gene-methylation Data}
\usage{
hf.snps.meth(B = 400, geno, meth, y, geno.pos, meth.pos, window.size,
  n.sample = nrow(geno), n.pair = 1000)
}
\arguments{
\item{B}{a numeric number specifying the number of bootstrapping times. Default is 400.}

\item{geno}{a data frame or matrix containing genotypes in the columns. Genotypes should be coded as (0, 1, 2) or (0, 1). SNP names should be stored as column names.}

\item{meth}{a data frame or matrix containing methylation data in the columns. Methylation data should be recoded as (0, 1, 2) or (0, 1). Names of CpG sites should be stored as column names.}

\item{y}{a numeric vector of 0 or 1.}

\item{geno.pos}{a data frame containing SNP names and positions in two columns.}

\item{meth.pos}{a data frame containing CpG names and positions in two columns.}

\item{window.size}{a numeric number specifying the size of genome distance. Interaction of the SNPs and CpG sites located within the size of genome distance will be evaluated exhaustively.}

\item{n.sample}{a numeric number specifying the number of samples to be included for estimating parameters. Default is the total number of samples.}

\item{n.pair}{a numeric value, the number of SNP-CpG pairs to use in bootstrapping. Default = min(P, 1000). P is the total number of pairs within the \code{window.size}.}
}
\value{
a set of \emph{h} and \emph{f} values indexed by \emph{k}, estimated automatically. Variable \emph{k} is the number of categorical combinations of a variable pair.
}
\description{
Estimate parameters (\emph{h} and \emph{f}) for \code{W-test}.
}
\examples{
data(SNP.pos)
data(CpG.pos)
data(genotype)
data(methylation)
data(phenotype2)

# Please note that parameter B is recommended to be greater than 400.
hf.pair <- hf.snps.meth(B = 80, geno = genotype, meth = methylation, y = phenotype2,
                        geno.pos = SNP.pos, meth.pos = CpG.pos, window.size = 1000)

}
\references{
Maggie Haitian Wang, Rui Sun, Junfeng Guo, Haoyi Weng, Jack Lee, Inchi Hu, Pak Sham and Benny C.Y. Zee (2016). A fast and powerful W-test for pairwise epistasis testing. Nucleic Acids Research. doi:10.1093/nar/gkw347.
}
\author{
Rui Sun, Maggie Haitian Wang
}
