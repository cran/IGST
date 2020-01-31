###########################################################
requireNamespace("e1071")
requireNamespace("BootMRMR")

IGST.pval.bootmrmrsvm<-function (x, y, re, Q, v)
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
  if (Q < 0 & Q > 1) {
    warning("Q is the quartile value of rank scores and must be within 0 and 1")
  }
  if (missing(Q)) {
    Q <- 0.5
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
    rankingCriteria <- v * w + (1-v) * (rsi/qsi)
    GeneRankedList <- sort(-rankingCriteria, index.return = TRUE)$ix
    rankvalue <- sort(GeneRankedList, index.return = TRUE)$ix
    rankscore <- (n1 + 1 - rankvalue)/(n1)
    M1[, j] <- as.vector(rankscore)
  }
  rankscore <- as.matrix(M1)
  mu <- Q
  R <- rankscore - mu
  sam <- nrow(R)
  pval.vec <- vector(mode = "numeric", length = nrow(rankscore))
  for (i in 1:sam) {
    z <- R[i, ]
    z <- z[z != 0]
    n11 <- length(z)
    r <- rank(abs(z))
    tplus <- sum(r[z > 0])
    etplus <- n11 * (n11 + 1)/4
    vtplus <- n11 * (n11 + 1) * (2 * n11 + 1)/24
    p.value = pnorm(tplus, etplus, sqrt(vtplus), lower.tail = FALSE)
    pval.vec[i] = p.value
  }

  class(pval.vec) <- "p values"
  return(pval.vec)
}
