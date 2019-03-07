library(bigchisqsum)

##--------------------------##
## code for Bausch's method ##
##--------------------------##


## Qa

lambda <- 200/(1:200)

## Bausch's method
N0 <- OptM(lambda[1], lambda[2], 5000, 0.5/50)

a <- chisqtwo(lambda[1], lambda[2], N = N0)

r <- vector("list", 100)
r[[1]] <- a
for (j in 1:99) {
  b <- chisqtwo(lambda[j * 2 + 1], lambda[j * 2 + 2], N = OptM(lambda[j * 2 + 1], 
                                                               lambda[j * 2 + 2], 5000, 0.1/50, memo = N0))
  a <- a %*% b
  print(c(sum(stcoef(a)), sum(sttail(a, 100))))
  r[[j + 1]] <- a <- simplify_kahan(a)
  print(j)
}

## compute relative error
q = c(1, 2, 5.3, 8, 9.5, 11, 12.5, 15) * 100
(dav = sapply(q, function(x) davies(x, lambda, rep(1, 200), acc = 1e-07, lim = 5e+07)$Qq))
bau = evaltail(r[[100]], q)
r = (bau - dav)/dav


## Qb

lambda <- 200/(1:20)

## Bausch's method
N0 <- OptM(lambda[1], lambda[2], 5000, 0.5/50)

a <- chisqtwo(lambda[1], lambda[2], N = N0)

r <- vector("list", 10)
r[[1]] <- a
for (j in 1:9) {
  b <- chisqtwo(lambda[j * 2 + 1], lambda[j * 2 + 2], N = OptM(lambda[j * 2 + 1], 
                                                               lambda[j * 2 + 2], 5000, 0.1/50, memo = N0))
  a <- a %*% b
  print(c(sum(stcoef(a)), sum(sttail(a, 100))))
  r[[j + 1]] <- a <- simplify_kahan(a)
  print(j)
}

## compute relative error
q2 = c(0.1, 0.2, 0.5, 2, 3.6, 5, 8.1, 10.5) * 100
(dav = sapply(q2, function(x) davies(x, lambda, rep(1, 20), acc = 1e-07, lim = 5e+07)$Qq))
bau = evaltail(r[[10]], q2)
r = (bau - dav)/dav