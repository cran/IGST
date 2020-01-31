\name{IGST.bootmrmrsvm.pval.cutoff}
\alias{IGST.bootmrmrsvm.pval.cutoff}
\title{
Gene selection using SVM and MRMR feature selection techniques with bootstrapping procedure
}
\description{
The informative gene set which have maximum relevance with target class or trait and minimum redundancy among genes based on statistical significance values computed from the combination of SVM and MRMR feature selection techniques with bootstrapping procedure.
}
\usage{
IGST.bootmrmrsvm.pval.cutoff (x, y, s, Q, v, re)
}

\arguments{
  \item{x}{ x is a n by p data frame of gene expression values where rows represent genes and columns represent samples. Each cell entry represents the expression level of a gene in a sample or subject (row names of x as gene names or gene ids).}
\item{y}{ y is a p by 1 numeric vector with entries 1 or -1 representing sample labels, where, 1\/-1 represents the sample label of subjects orsamples for stress/control condition(for two class problems).}
\item{s}{s is a scalar representing the size of the informative gene set to be obtained.}
\item{Q}{Q is a scalar representing the quartile value of the rank scores of genes (lies within 		1\/N to 1), usually the second quartile, i.e. 0.5 or third quartile i.e. 0.75 may be taken.}
\item{v}{v is a scalar representing the weightage of a method and must be within 0 and 1. }
\item{re}{re is a scalar representing the number of bootstrap generated, re must be sufficiently large (i.e. number of times bootstrap samples are generated.}
}

\value{The function returns a list of the genes\/informative gene set which are highly relevant to the particular trait or condition under investigation and minimal redundant among themselves without any spurious association among the genes.}

\references{
Das, S., Rai, A., Mishra, D. C., & Rai, S. N. (2018). Statistical approach for selection of biologically informative genes. Gene, 655, 71-83.

Guyon, I., Weston, J., Barnhill, S., & Vapnik, V. (2002). Gene selection for cancer classification using support vector machines. Machine Learning, 46(1-3), 389-422.

Wang, J., Chen, L., Wang, Y., Zhang, J., Liang, Y., & Xu, D. (2013). A computational systems biology study for understanding salt tolerance mechanism in rice. PLoS One, 8(6), e64929.

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
Q<-0.5
v<-0.5
re<-10
IGST.bootmrmrsvm.pval.cutoff (x, y, s, Q, v, re)
}
