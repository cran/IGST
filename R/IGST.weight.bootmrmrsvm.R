#################################################################
requireNamespace("e1071")
requireNamespace("BootMRMR")

IGST.weight.bootmrmrsvm<-function (x, y, re, v)
{

  this.call = match.call()
  if ((!class(x) == "data.frame")) {
    warning("x must be a data frame and rows as gene names")
  }
  if ((!class(y) == "numeric")) {
    warning("y must be a vector of 1/-1's for two class problems")
  }
  if (!length(y) == ncol(x)) {
    warning("Number of samples in x must have same number of sample labels in y")
  }
  if (re < 0 & re <= 50) {
    warning("s must be numeric and sufficiently large")
  }

  if (v < 0 & v > 1) {
    warning("v is the tradeoff value between svm and mrmr and must be within 0 and 1 ")
  }



  cls <- as.numeric(y)
  genes <- rownames(x)
  g <- as.matrix(x)
  n1 <- nrow(g)
  M <- ncol(x)
  GeneRankedList <- vector(length = n1)
  M1 <- matrix(0, n1, re)
  for (j in 1:re) {
    samp <- sample(M, M, replace = TRUE)
    x1 <- g[, samp]
    y1 <- cls[samp]
    qsi <- as.vector((apply(abs(cor(t(x1), method = "pearson",
                                    use = "p") - diag(n1)), 1, sum))/(n1 - 1))
    idx <- which(y1 == 1)
    idy <- which(y1 == -1)
    B = vector(mode = "numeric", n1)
    for (i in 1:nrow(x1)) {
      f.mes <- (((mean(x1[i, idx]) - mean(x1[i, ]))^2) +
                  ((mean(x1[i, idy]) - mean(x1[i, ]))^2))/(var(x1[i,
                                                                  idx]) + var(x1[i, idy]))
      B[i] <- f.mes
    }
    svmModeli = svm(t(x1), as.matrix(y1), cost = 10,
                    cachesize = 500, scale = FALSE, type = "C-classification",
                    kernel = "linear")
    w = abs(as.vector(t(svmModeli$coefs) %*% svmModeli$SV))
    rsi <- abs(B)
    Weight <- v * w + (1-v) * (rsi/qsi)

    GeneRankedList <- sort(-Weight, index.return = TRUE)$ix

    rankvalue <- sort(GeneRankedList, index.return = TRUE)$ix
    rankscore <- (n1 + 1 - rankvalue)/(n1)
    M1[, j] <- as.vector(rankscore)
  }
  #rankingCriteria=Weight
  class(Weight) <- "Weight values"
  return(Weight)
}
