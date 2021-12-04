#' Modified Kruskal Wallis test (ZIKW) for zero-inflated data
#'
#' @param x a numeric vector of data values with groups defined by the group parameter, or a list of numeric data vectors with each vector in the list as one group.
#' @param group a vector giving the groups for the corresponding elements of \code{x}. Ignored if \code{x} is a list.
#' @param perm use permutations to calculate the p-value
#' @param perm.n number of permutations used for p-value calculation
#' @return modified Kruskal Wallis test statistic and p-value
#' @export
#' @import stats utils graphics methods
#' @examples
#' ## x is a list
#' x <- list(group1 = c(rep(0, 5),rlnorm(20, meanlog = 0, sdlog = 1)),
#'      group2=c(rep(0, 10),rlnorm(20, meanlog = 1, sdlog = 1)),
#'      group3=c(rep(0, 15),rlnorm(20, meanlog = 2, sdlog = 1)))
#' zikw(x, perm = FALSE)
#' ## x is a vector
#' x <- c(c(rep(0, 5),rlnorm(20, meanlog = 0, sdlog = 1)),
#'     c(rep(0, 10),rlnorm(20, meanlog = 1, sdlog = 1)),
#'     c(rep(0, 15),rlnorm(20, meanlog = 2, sdlog = 1)))
#' group <- c(rep('group1', 25),rep('group2', 30),rep('group3', 35))
#' zikw(x, group, perm = FALSE)
#' ## use permutations to calculate the pvalue
#' \dontrun{zikw(x, group, perm = TRUE)}
#' @references Wanjie Wang, Eric Z. Chen and Hongzhe Li (2021). Rank-based tests for compositional distributions with a clump of zeros. Submitted.

zikw <- function(x, group, perm = FALSE, perm.n = 10000) {

  ## check if x is a list
  if (class(x) == "numeric" || class(x) == "data.frame") {
    x <- vect2list(x, group)
  }
  
  ## calclulate the zikw statistic
  w <- calculate_zikw_statistic(x)
  
  ## calculate the pvalue
  if (perm) {
    permu.w <- rep(0, perm.n)
    perm.group <- rep(c(seq_len(x)), sapply(x, length))
    xvec <- unlist(x)
    ## need to permute a list
    for (i in seq_len(perm.n)) {
      set.seed(i)
      xvec.perm <- sample(xvec, length(xvec))
      x.list.perm <- vect2list(xvec.perm, group)
      permu.w[i] <- calculate_zikw_statistic(x.list.perm)
    }
    pw <- sum(abs(w) < abs(permu.w)) / perm.n
  } else {
    K <- length(x)
    pw <- pchisq(w, K - 1, lower.tail = F)
  }
  
  return(list(p.value = pw, statistics = w))
}

calculate_zikw_statistic <- function(x) {
  
  ## number of groups
  K <- length(x)
  
  ## total observations in each group
  N <- rep(0, K)
  
  ## number of non-zero observations in each group
  n <- rep(0, K)
  xvec <- numeric(0)
  
  ## count total and non-zero observations in each group
  for (i in seq_len(K)) {
    N[i] <- length(x[[i]])
    n[i] <- sum(x[[i]] != 0)
    xvec <- c(xvec, x[[i]])
  }
  
  ## non-zero proportion
  prop <- n / N
  pmax <- max(prop)
  
  ## keep only round(pmax * N) observations in each group
  Ntrun <- round(pmax * N)
  
  ## truncate zeros in each group
  Xtrun.vec <- numeric(0)
  for (i in seq_len(K)) {
    data <- x[[i]]
    Xtrun.vec <-
      c(Xtrun.vec, data[data != 0], rep(0, Ntrun[i] - n[i]))
  }
  rankdata <- sum(Ntrun) + 1 - rank(Xtrun.vec)
  r <- sum(rankdata[seq_len(Ntrun[i])])
  
  for (i in seq_len(K)[-1]) {
    index <- seq_len(Ntrun[i]) + sum(Ntrun[seq_len(i - 1)])
    r <- c(r, sum(rankdata[index]))
  }
  s <- r - Ntrun * (sum(Ntrun) + 1) / 2
  u <- 0L
  
  for (i in seq_len(K - 1))
    u <- c(u, N[i + 1] * sum(s[1:i]) - sum(N[1:i]) * s[i + 1])
  u <- u / sum(N)^2
  thetam <- mean(prop)
  simun <- matrix(0, nrow = 5000, ncol = K)
  simup <- simun
  for (ss in seq_len(K)) {
    simun[,ss] <- rbinom(5000, N[ss], thetam)
    simup[,ss] <- simun[,ss] / N[ss]
  }
  simupmax <- apply(simup, 1, max)
  varsimu <- numeric(K - 1)
  
  varsimu[1] <-
    N[2]^2 * mean(simupmax^2 * (simup[,1] - simup[,2])^2) * N[1]^2
  varu2 <- N[2] * N[1] * (N[1] + N[2])
  
  for (ss in seq_len(K - 1)[-1]) {
    varsimu[ss] <-
      N[ss + 1]^2 * mean(simupmax^2 * (apply(simun[,seq_len(ss)], 1, sum) - 
                                         simup[,ss + 1] * sum(N[seq_len(ss)])) ^ 2)
    varu2 <-
      c(varu2, N[ss + 1] * sum(N[seq_len(ss)]) * sum(N[seq_len(ss + 1)]))
  }
  varsimu <- varsimu / (sum(N))^2 / 4
  
  varu2 <-
    varu2 * thetam^2 * (thetam + 1 / sum(N)) / 12 / (sum(N))^2
  varu <- varsimu + varu2
  
  ## modified Kruskal Wallis test statistic
  w <- sum(u^2 / varu)
  return(w)
}

vect2list <- function(x, group){
  s <- unique(group)
  x.list <- list()
  for (n in s) {
    x.list[[n]] <- x[group == n]
  }
  return(x.list)
}
