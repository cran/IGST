\name{IGST.bootmrmrsvm.weight.cutoff}
\alias{IGST.bootmrmrsvm.weight.cutoff}
\title{
Identification of informative gene set based on weights obtained from SVM and MRMR feature selection technique with bootstrapping procedure
}
\description{
The function enables to find set of informative genes based on weights which are obtained by maximizing the relevancy of genes with classes or condition or trait and minimizing the redundancy among genes from the combination of SVM and MRMR feature selection techniques with bootstrapping procedure.
}
\usage{
IGST.bootmrmrsvm.weight.cutoff (x, y, s, v, re)
}

\arguments{
  \item{x}{ x is a n by p data frame of gene expression values where rows represent genes and columns represent samples. Each cell entry represents the expression level of a gene in a sample or subject (row names of x as gene names or gene ids).}
\item{y}{ y is a p by 1 numeric vector with entries 1 or -1 representing sample labels, where, 1 or -1 represents the sample label of subjects or samples for stress or control condition(for two class problems).}
\item{s}{s is a scalar representing the size of the informative gene set to be obtained.}
\item{v}{v is a scalar representing the weightage of a method and must be within 0 and 1. }
\item{re}{re is a scalar representing the number of bootstrap generated, re must be sufficiently large (i.e. number of times bootstrap samples are generated.}
}

\value{The function returns a set of genes, which are highly informative to the trait or condition under consideration based on weights given by the combination of SVM and MRMR feature selection techniques with bootstrapping procedure. }
\references{
Ding, C., & Peng, H. (2005). Minimum redundancy feature selection from microarray gene expression data. Journal of Bioinformatics and Computational Biology, 3(02), 185-205.

Mishra DC, Kumar S, Lal SB, Saha A, Chaturvedi KK, Budhlakoti N, et al.( 2018) TAGPT: A Web Server for Prediction of Trait Associated Genes using Gene Expression Data. Annals of Genetics and Genetic Disorder. 1(1): 1003.
}
\author{
Nitesh Kumar Sharma, Dwijesh Chandra Mishra, Neeraj Budhlakoti and Md. Samir Farooqi
}

\examples{
#################################
library(IGST)
data(rice_cold)
x<-rice_cold[-1,]
y<-rice_cold[1,]
y<-as.matrix(y)
y<-as.vector(y)
s<-10
#Q<-0.5
v<-0.5
re<-10
IGST.bootmrmrsvm.weight.cutoff (x, y, s, v, re)
}
