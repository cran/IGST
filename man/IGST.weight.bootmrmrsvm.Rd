\name{IGST.weight.bootmrmrsvm}
\alias{IGST.weight.bootmrmrsvm}
\title{
Computation of weights for informative genes or gene set selection using SVM and MRMR feature selection technique with bootstrapping procedure
}
\description{
The function computes the weights associated with each genes for a given dataset using SVM and MRMR feature selection technique with bootstrapping procedure.
}
\usage{
IGST.weight.bootmrmrsvm (x, y, re, v)
}
\arguments{
  \item{x}{ x is a n by p data frame of gene expression values where rows represent genes and columns represent samples. Each cell entry represents the expression level of a gene in a sample or subject (row names of x as gene names or gene ids).}
\item{y}{ y is a p by 1 numeric vector with entries 1 or -1 representing sample labels, where, 1\/-1 represents the sample label of subjects or samples for stress or control condition(for two class problems).}
\item{v}{v is a scalar representing the weightage of a method and must be within 0 and 1. }
\item{re}{re is a scalar representing the number of bootstrap generated, re must be sufficiently large (i.e. number of times bootstrap samples are generated.}
}

\value{The function returns a vector of weights associated with each genes computed from SVM and MRMR feature selection technique with bootstrapping procedure for a given dataset. }
\references{
Wang, J., Chen, L., Wang, Y., Zhang, J., Liang, Y., & Xu, D. (2013). A computational systems biology study for understanding salt tolerance mechanism in rice. PLoS One, 8(6), e64929.

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
#s<-10
#Q<-0.5
v<-0.5
re<-10
IGST.weight.bootmrmrsvm (x, y, re, v)
}
