\name{IGST.pval.bootmrmrsvm}
\alias{IGST.pval.bootmrmrsvm}
\title{
Computation of statistical significance values for genes using SVM and MRMR feature selection technique with bootstrapping procedure
}
\description{
The function computes the statistical significance values for the genes from the non-parametric test "H0: i-th gene is not informative against H1: i-th gene is informative" for selection of informative genes using SVM and MRMR feature selection technique with bootstrapping procedure.
}
\usage{
IGST.pval.bootmrmrsvm(x, y, re, Q, v)
}

\arguments{
  \item{x}{ x is a n by p data frame of gene expression values where rows represent genes and columns represent samples. Each cell entry represents the expression level of a gene in a sample or subject (row names of x as gene names or gene ids).}
\item{y}{ y is a p by 1 numeric vector with entries 1 or -1 representing sample labels, where, 1 or -1 represents the sample label of subjects or samples for stress or control condition(for two class problems).}
\item{Q}{Q is a scalar representing the quartile value of the rank scores of genes (lies within 1\/N to 1), usually the second quartile, i.e. 0.5 or third quartile i.e. 0.75 may be taken.}
\item{v}{v is a scalar representing the weightage of a method and must be within 0 and 1. }
\item{re}{re is a scalar representing the number of bootstrap generated, re must be sufficiently large (i.e. number of times bootstrap samples are generated.}
}

\value{
The function returns a vector of p-values for all the genes from the given statistical test in the dataset using SVM and MRMR feature selection technique with bootstrapping procedure.
}
\references{
Das, S., Rai, A., Mishra, D. C., & Rai, S. N. (2018). Statistical approach for selection of biologically informative genes. Gene, 655, 71-83.
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
Q<-0.5
v<-0.5
re<-10
IGST.pval.bootmrmrsvm (x, y, re, Q, v)
}


