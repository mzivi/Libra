require(reshape2)
require(ggplot2)

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
  
  # state tracks the probabilities for each n: state[n+1] = P[k|n]
  state <- array(NA, n+1)
  state[1:min(k+1, n+1)] <- 1 # don't need to calculate for n <= k
  
  CalcPLongestRunLessEqualImpl1(p, k, n, state)
}

CalcPLongestRunLessEqualImpl1 <- function(p, k, n, state) {

  if (n <= k) # we have alerady dealt with these during initialization
    return (state)
  
  pkn <- state[n+1]
  if (is.na(pkn)) { # it's a first timer!
    
    # this will update state up to n-1
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
  s <- c(1, rep(0, k+1)) # initial state, k = 0
  result[1] <- 1 # n = 0
  if (n > 0) {
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

kmax <- 10
n <- 10
p <- 0.5

Prob <- matrix(NA, kmax + 1, n + 1)
for (k in 0:kmax) {
  Prob[k+1,] <- 1-CalcPLongestRunLessEqual2(p, k, n)
}
rownames(Prob) <- 0:(nrow(Prob)-1)
colnames(Prob) <- 0:(ncol(Prob)-1)

Prob <- melt(Prob)
colnames(Prob) <- c("k", "n", "P")
head(Prob[Prob$k==0, ])
head(Prob[Prob$k==1, ])
head(Prob[Prob$n==2, ])
head(Prob[Prob$n==3, ])

ggplot(data = Prob, aes(x=n, y=P, group=k, colour=factor(k))) +
  theme_bw() + ggtitle("P[ln > k| n") + xlab("n") + geom_line() + geom_point()


# let's tackle P[ln = k|n]
# just use approach 2 above, but when there has been a run of length k
# if next trial is tail jump into a second set of states which counts
# again the run length to make sure there is not going to be a longer one.
# state s = (s_a, s_b, T)
# s_a[k+1] = prob length of last run equals k conditional on there hasn't
# been a run of length k or more
# s_b[k=1] = prob length of last run equals k conditional on there has
# been a run of length k but not more.
# T = there has been a run of length more than k
CalcPLongestRunEqualK = function(p, k, n) {
  if (n < 0)
    stop("n < 0")
  if (k < 0)
    stop("k < 0")
  if (k == 0)
    return (rep(1 - p, n + 1)^(0:n))
  
  result <- rep(NA, n + 1)
  Aa <- rbind(rep(1 - p, k), cbind(p * diag(k - 1), rep(0, k - 1)))
  A.b.a <- matrix(0, k, k + 2) # from b no way to get back to a
  A.a.b <- matrix(0, k + 2, k)
  A.a.b[k + 1, k] <- p # only way from a -> b, first run reaches length k
  Ab <- cbind(rbind(rep(1 - p, k+1), p * diag(k+1)), c(rep(0, k+1), 1))
  A <- rbind(cbind(Aa, A.b.a), cbind(A.a.b, Ab))

  s <- c(1, rep(0, 2 * k + 1)) # initial state, k = 0
  result[1] <- ifelse(k == 0, 1, 0) # n = 0
  if (n > 0) {
    for (j in 1:n) {
      s <- A %*% s
      result[j+1] <- sum(s[(k + 1):(nrow(A) - 1)])
    }
  }
  result
}

CalcPLongestRunEqualKNaive = function(p, k, n) {
  if (n < 0)
    stop("n < 0")
  if (k < 0)
    stop("k < 0")
  if (k == 0)
    return (rep(1 - p, n + 1)^(0:n))
  if (n == 0)
    return (c(0))
  
  result <- rep(NA, n + 1)
  result[1] <- 0
  # init for n == 1
  for (i in 1:n) {
    if (i == 1) {
      all.sequences <- matrix(c(0, 1), 2, 1)
    } else {
      all.sequences <- rbind(cbind(rep(0, nrow(all.sequences)), all.sequences),
                             cbind(rep(1, nrow(all.sequences)), all.sequences))
    }
    n.k.max <- 0
    for (j in 1:nrow(all.sequences)) {
      max.last.length <- 0
      max.length <- 0
      for (l in 1:i) {
        if (all.sequences[j, l] == 1) {
          max.last.length = max.last.length + 1
          if (max.last.length > max.length)
            max.length <- max.last.length
        } else {
          max.last.length <- 0
        }
      }
      if (max.length == k)
        n.k.max <- n.k.max + 1
    }
    result[i + 1] <- n.k.max / nrow(all.sequences)
  }

  result
}


kmax <- 10
n <- 100
p <- 0.5
Prob <- matrix(NA, kmax + 1, n + 1)
for (k in 0:kmax) {
  Prob[k+1,] <- CalcPLongestRunEqualK(p, k, n)
}
rownames(Prob) <- 0:(nrow(Prob)-1)
colnames(Prob) <- 0:(ncol(Prob)-1)


# let's check a few...
kmax <- 10
n <- 20
p <- 0.5
ProbNaive <- matrix(NA, kmax + 1, n + 1)
for (k in 0:kmax) {
  ProbNaive[k+1,] <- CalcPLongestRunEqualKNaive(p, k, n)
}
rownames(ProbNaive) <- 0:(nrow(ProbNaive)-1)
colnames(ProbNaive) <- 0:(ncol(ProbNaive)-1)

# they should be exactly the same
sum(ProbNaive - Prob[1:nrow(ProbNaive), 1:ncol(ProbNaive)] )

# but the Markov Chain is clearly much faster
system.time(result1 <- CalcPLongestRunEqualK(0.5, 1, 18))
system.time(result2 <- CalcPLongestRunEqualKNaive(0.5, 1, 18))

Prob <- melt(Prob)
colnames(Prob) <- c("k", "n", "P")
head(Prob[Prob$k==0, ])
head(Prob[Prob$k==1, ])
head(Prob[Prob$n==2, ])
head(Prob[Prob$n==3, ])

ggplot(data = Prob, aes(x=n, y=P, group=k, colour=factor(k))) +
  theme_bw() + ggtitle("P[ln = k| n") + xlab("n") + geom_line() + geom_point()
