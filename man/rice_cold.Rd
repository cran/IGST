\name{rice_cold}
\alias{rice_cold}
\docType{data}
\title{
A gene expression dataset of rice under cold stress
}
\description{
This data has gene expression values of 250 genes over 36 samples or subjects for a cold vs. control study in rice. These 36 samples belong to either of cold stress or control condition (two class problem). This gene expression data is balanced type as the first 18 samples are under cold stress and the later 18 samples are under control condition. The first row of the data contains the samples or subjects labels with entries are 1 and -1, where the label ’1’ and ’-1’ represent samples generated under cold stress and control condition respectively.
}
\usage{data("rice_cold")}
\format{
  A data frame with 250 rows as genes with 36 columns as samples or subjects. Each column (sample) represent the gene expression values of genes. Each column as microarray samples with labels -1 or 1 represents control or cold stress respectively.
}
\details{
The data is created by taking 250 genes from the large number of genes from NCBI GEO database. The rows are the genes and columns are the samples or subjects. The first half of the samples or subjects are generated under cold stress condition and other half under control condition. The first row of the data contains the samples/subjects labels with entries are 1 and -1, where the label ’1’ and ’-1’ represent samples generated under cold stress and control condition respectively.
}
\source{
Gene Expression Omnibus: NCBI gene expression and hybridization array data repository.ncbi.nlm.nih.gov/geo/.
}

\examples{
####################################
library(IGST)
data(rice_cold)

}

