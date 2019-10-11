

## Date create: March 07, 2019

##------------------##
## code for Table 2 ##
##------------------##

library(CompQuadForm)
library(bigQF)
library(survey)

## Currently, if you use pQF function in bigQF package and choose Davies's method
## to combine the leading eigenvalues, you cannot get same result as Table 2,
## because the accuracy is set to be $10^-9$ in the function. In table 2, we change
## it to $10^{-16}$.  You may get a larger error when p-value is greater than
## $10^{-9}$ if you run the code directly without changing accuracy. But you can get 
## the same result as Table 2 if you modify the accuracy to be $10^{-16}$. This can 
## be done by set it as a global argument.

## set global accuracy for pQF funtion (the defalut value is 10^{-6})
acc = 1e-16

set.seed(201902)

## Q1

## read the matrix and non-zero eigenvalues
X = readRDS("s500.rds")
eigen = readRDS("eigenVal500.rds")
q1 = seq(12, length.out = 4, by = 3.6) * 10^3

## compute the relative error
(dav1 = sapply(q1, function(x) davies(x, eigen, rep(1, 305), acc = 1e-16, lim = 5e+07)$Qq))
(fast1 = pQF(q1, X, neig = 50, convolution.method = "integration"))
(error1 = (fast1 - dav1)/dav1)


## Q2

## read the matrix and non-zero eigenvalues
X = readRDS("s1000.rds")
eigen = readRDS("eigenVal1000.rds")
q2 = seq(4, length.out = 4, by = 1.4) * 10^4

## compute the relative error
(dav2 = sapply(q2, function(x) davies(x, eigen, rep(1, 637), acc = 1e-16, lim = 5e+07)$Qq))
(fast2 = pQF(q2, X, neig = 50, convolution.method = "integration"))
(error2 = (fast2 - dav2)/dav2)



## Q3

## read the matrix and non-zero eigenvalues
X = readRDS("s2000.rds")
eigen = readRDS("eigenVal2000.rds")
q3 = seq(11, length.out = 4, by = 4) * 10^4

## compute the relative error
(dav3 = sapply(q3, function(x) davies(x, eigen, rep(1, 1063), acc = 1e-16, lim = 5e+07)$Qq))
(fast3 = pQF(q3, X, neig = 100, convolution.method = "integration", tr2.sample.size = 2000))
(error3 = (fast3 - dav3)/dav3)



## Q4

## read the matrix and non-zero eigenvalues
X = readRDS("s7000.rds")
eigen = readRDS("eigenVal7000.rds")
q4 = seq(1.2, length.out = 4, by = 0.5) * 10^6

## compute the relative error
(dav4 = sapply(q4, function(x) davies(x, eigen, rep(1, 3985), acc = 1e-16, lim = 5e+07)$Qq))
(fast4 = pQF(q4, X, neig = 100, convolution.method = "integration"))
(error4 = (fast4 - dav4)/dav4)


## Q5

## read the matrix and non-zero eigenvalues
X = readRDS("s9000.rds")
eigen = readRDS("eigenVal9000.rds")
q5 = seq(2, length.out = 4, by = 0.5) * 10^6

## compute the relative error
(dav5 = sapply(q5, function(x) davies(x, eigen, rep(1, 4984), acc = 1e-16, lim = 5e+07)$Qq))
(fast5 = pQF(q5, X, neig = 200, convolution.method = "integration"))
(error5 = (fast5 - dav5)/dav5)


## Q6

## read the matrix and non-zero eigenvalues
eigen = readRDS("eigenVal20000.rds")
q6 = seq(9, length.out = 4, by = 3) * 10^6

## compute the relative error
(dav6 = sapply(q6, function(x) davies(x, eigen, rep(1, 11259), acc = 1e-16, lim = 5e+07)$Qq))
(fast6 = pQF(q6, diag(eigen), neig = 200, convolution.method = "integration"))
(error6 = (fast6 - dav6)/dav6)



##------------------##
## code for Table 3 ##
##------------------##

# Q1

## read the matrix and non-zero eigenvalues
set.seed(201902)
X = readRDS("s500.rds")
eigen = readRDS("eigenVal500.rds")
q1 = seq(28, 58, by = 10) * 10^3

## compute the relative error
(saddle1 = pchisqsum(q1, lower.tail = FALSE, df = rep(1, 305), a = eigen, method = "sad"))
(fasts1 = pQF(q1, X, neig = 50))
(error1 = (fasts1 - saddle1)/saddle1)

## Q2

## read the matrix and non-zero eigenvalues
X = readRDS("s1000.rds")
eigen = readRDS("eigenVal1000.rds")
q2 = seq(8, length.out = 4, by = 2) * 10^4

## compute the relative error
(saddle2 = pchisqsum(q2, lower.tail = FALSE, df = rep(1, 637), a = eigen, method = "sad"))
(fasts2 = pQF(q2, X, neig = 50))
(error2 = (fasts2 - saddle2)/saddle2)


## Q3

## read the matrix and non-zero eigenvalues
X = readRDS("s2000.rds")
eigen = readRDS("eigenVal2000.rds")
q3 = seq(2, length.out = 4, by = 0.5) * 10^5

## compute the relative error
(saddle3 = pchisqsum(q3, lower.tail = FALSE, df = rep(1, 1063), a = eigen, method = "sad"))
(fasts3 = pQF(q3, X, neig = 100, tr2.sample.size = 2000))
(error3 = (fasts3 - saddle3)/saddle3)

## Q4

## read the matrix and non-zero eigenvalues
X = readRDS("s7000.rds")
eigen = readRDS("eigenVal7000.rds")
q4 = seq(3, length.out = 4, by = 0.5) * 10^6

## compute the relative error
(saddle4 = pchisqsum(q4, lower.tail = FALSE, df = rep(1, 3985), a = eigen, method = "sad"))
(fasts4 = pQF(q4, X, neig = 100))
(error4 = (fasts4 - saddle4)/saddle4)


## Q5

## read the matrix and non-zero eigenvalues
eigen = readRDS("eigenVal9000.rds")
q5 = seq(4, length.out = 4, by = 0.5) * 10^6

## compute the relative error
(saddle5 = pchisqsum(q5, lower.tail = FALSE, df = rep(1, 4984), a = eigen, method = "sad"))
(fasts5 = pQF(q5, diag(eigen), neig = 200))
(error5 = (fasts5 - saddle5)/saddle5)


# Q6

## read the matrix and non-zero eigenvalues
eigen = readRDS("eigenVal20000.rds")
q6 = seq(2, length.out = 4, by = 0.5) * 10^7

## compute the relative error
(saddle6 = pchisqsum(q6, lower.tail = FALSE, df = rep(1, 11259), a = eigen, method = "sad"))
(fasts6 = pQF(q6, diag(eigen), neig = 200))
(error6 = (fasts6 - saddle6)/saddle6)



##------------------##
## code for Table 4 ##
##------------------##

## case A
eigen = readRDS("eigenVal9000.rds")
q = 2300000
system.time(result <- davies(q, eigen, rep(1, 4984), acc = 1e-13, lim = 5e+07)$Qq)

## case B
eigen1 = eigen
eigen1[1] = eigen1[1] * 10
q = 1.3e+07
system.time(result <- davies(q, eigen1, rep(1, 4984), acc = 1e-13, lim = 5e+07)$Qq)

## case C
eigen2 = eigen1
eigen2[1] = eigen2[1] * 10
q = 1.2e+08
system.time(result <- davies(q, eigen2, rep(1, 4984), acc = 1e-13, lim = 5e+07)$Qq)

## case D
eigen3 = eigen2
eigen3[1] = eigen3[1] * 10
q = 1.2e+09
system.time(result <- davies(q, eigen3, rep(1, 4984), acc = 1e-13, lim = 5e+07)$Qq)

## case E
eigen4 = eigen3
eigen4[1] = eigen4[1] * 10
q = 1.2e+10
system.time(result <- davies(q, eigen4, rep(1, 4984), acc = 1e-13, lim = 5e+07)$Qq)


##------------------##
## code for Table 5 ##
##------------------##


## compare computational time of SVD and SSVD

## Q0
system.time(ssvd(SKAT.example$Z, n = 50))
system.time(svd(SKAT.example$Z)$d)

## Q1
X = readRDS("s500.rds")
system.time(ssvd(X, n = 100))
system.time(svd(X)$d)

## Q2
X = readRDS("s1000.rds")
system.time(ssvd(X, n = 100))
system.time(svd(X)$d)

## Q3
X = readRDS("s2000.rds")
system.time(ssvd(X, n = 100))
system.time(svd(X)$d)

## Q4
X = readRDS("s7000.rds")
system.time(ssvd(X, n = 100))
system.time(svd(X)$d)
