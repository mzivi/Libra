# Estimation of P[ln > k|n] where ln is the length of the longest run
# of successes in n Bernoulli trials.
#
# we start by calculating P[ln <= k|n] := P[k|n]

# Approach 1: write an equation for P[k|n] by iteratively
# contitioning from the first throw and on:
# P[k|n] = p P[k|n, H] + (1-p) P[k|n-1] =
# ...    = p^k (1-p) P[k|n-1-k] + ... + (1-p) P[k|n-1]
# solve this by dynamic programming:

CalcPLongestRunLessEqual1 <- function(p, k, n) {
  
  if (n < 0)
    stop("n < 0")
  if (k < 0)
    stop("k < 0")
  
  state <- array(NA, n+1)
  state[1:min(k+1, n+1)] <- 1
  CalcPLongestRunLessEqualImpl1(p, k, n, state)
}

CalcPLongestRunLessEqualImpl1 <- function(p, k, n, state) {

  if (n <= k)
    return (state)
  
  pkn <- state[n+1]
  if (is.na(pkn)) {
    state <- CalcPLongestRunLessEqualImpl1(p, k, n - 1, state)
    pkn = 0
    for (j in 0:min(k, n-1)) {
      pkn = pkn + p^j * state[n - j]
    }
    pkn <- (1-p) * pkn
    state[n+1] <- pkn
  }
  state
}

# Approach 2: use a Markov Chain for the evolution of state variable
# s = 0, 1, 2, ..., k length of the latest run + a final absorbing
# state (s[k+2]) which embed when there has been a run longer than k.
# The transition matrix A has the form:
#
# 1-p 1-p 1-p ... 1-p  0
#   p   0   0 ...   0  0
#   0   p   0 ...   0  0
#   0   0   p ...   0  0
#   ...
#   0   0   0       p  1
#
# given the initial state s = (1, 0, 0, ..., 0), the state probs
# after n throws is A^n, therefore P[k|n] being the probability
# not to be in the final absorbing state (run longer than k), we have
# P[k|n] = (1 - A^n s)[k+2]

CalcPLongestRunLessEqual2 = function(p, k, n) {
  if (n < 0)
    stop("n < 0")
  if (k < 0)
    stop("k < 0")
  
  result <- rep(NA, n+1)
  A <- cbind(rbind(rep(1-p, k+1), p * diag(k+1)), c(rep(0, k+1), 1))
  s <- c(1, rep(0, k+1))
  result[1] <- 1 # n = 0
  if (n>0) {
    for (j in 1:n) {
      s <- A %*% s
      result[j+1] <- 1-s[k+2]
    }
  }
  result
}

# approach 2 is slightly faster:
system.time(result1 <- CalcPLongestRunLessEqual1(0.5, 1, 500))
system.time(result2 <- CalcPLongestRunLessEqual2(0.5, 1, 500))

# also, approach one runs out of stacks!
system.time(result1long <- tryCatch(
  CalcPLongestRunLessEqual1(0.5, 100, 100000),
  error = function(c) message("Too many recursive calls.")))
# while the Markov Chain can go on a lot
system.time(result2long <- try(CalcPLongestRunLessEqual2(0.5, 100, 100000)))

# OK, they give the same answer:
sum(result1-result2)

# We now want to calculate P[ln>k|n], which is just 1 - P[k|n]
# as defined above.

Prob <- matrix(NA, 11, 51)
for (k in 0:10) {
  Prob[k+1,] <- 1-CalcPLongestRunLessEqual2(0.5, k, 50)
}
rownames(Prob) <- 0:(nrow(Prob)-1)
colnames(Prob) <- 0:(ncol(Prob)-1)

if (!require(reshape2)) {
  head(Prob)
  stop("Please install reshape2 package in order to continue.")
}

Prob <- melt(Prob)
colnames(Prob) <- c("k", "n", "P")
head(Prob[Prob$k==0, ])
head(Prob[Prob$k==1, ])
head(Prob[Prob$n==2, ])
head(Prob[Prob$n==3, ])

if (require(ggplot2)) {
  p = ggplot(data = Prob, aes(x=n, y=P, group=k, colour=factor(k))) +
    theme_bw() + ggtitle("P[ln > k| n") + xlab("n") + geom_line() + geom_point()
  plot(p)
}
