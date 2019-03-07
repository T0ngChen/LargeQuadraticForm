
library(CompQuadForm)
library(survey)

##-------------------##
## code for Figure 1 ##
##-------------------##


## Q1 
## read non-zero eigenvalues
eigen = readRDS("eigenVal500.rds")
q = (c(seq(6.6, by = 1.9, length.out = 3), 12, 14, 17.5, 21, 24.5, 28.5) * 10^3)

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 305), acc = 1e-16, lim = 5e+07)$Qq))

## Farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 305), eps = 1e-16, maxit = 2.14e+09)$Qq))

## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 305))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 305), a = eigen, method = "satterthwaite"))

## saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 305), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## plot
pdf("Q1.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of error ratio", 
     main = expression(Q[1]), xaxt = "n", cex.main = 1.4, cex.lab = 1.4, yaxt = "n")
lines(q, saddle_relative, cex = 1.2, lty = 4)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.4 %*% 10^-1), expression(1.3 %*% 10^-3), 
                                                  expression(1.2 %*% 10^-5), expression(1.3 %*% 10^-7), expression(1.5 %*% 
                                                                                                                                                                     10^-9), expression(1.8 %*% 10^-11), expression(1.1 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite", 
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()


## Q2
## read non-zero eigenvalues
eigen = readRDS("eigenVal1000.rds")
q = (c(seq(2.3, by = 0.55, length.out = 5), 5.6, 6.7, 7.8, 9) * 10^4)

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 637), acc = 1e-16, lim = 5e+07)$Qq))

## Farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 637), eps = 1e-16, maxit = 2.14e+09)$Qq))

## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 637))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 637), a = eigen, method = "satterthwaite"))

## saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 637), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## plot
pdf("Q2.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of error ratio", 
     main = expression(Q[2]), xaxt = "n", cex.main = 1.4, cex.lab = 1.4, yaxt = "n")
lines(q, saddle_relative, cex = 1.2, lty = 4)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.4 %*% 10^-1), expression(1.6 %*% 10^-3), 
                                                  expression(1.4 %*% 10^-5), expression(1.4 %*% 10^-7), expression(1.5 %*% 
                                                                                                                                                                     10^-9), expression(1.7 %*% 10^-11), expression(1.3 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite", 
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()




## Q3
## read non-zero eigenvalues
eigen = readRDS("eigenVal2000.rds")
q = (c(seq(6.7, by = 1.5, length.out = 5), 15.5, 18.5, 21.5, 24.5) * 10^4)

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 1063), acc = 1e-16, lim = 5e+07)$Qq))

## Farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 1063), eps = 1e-16, maxit = 2.14e+09)$Qq))

## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 1063))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 1063), a = eigen, method = "satterthwaite"))

## saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 1063), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## Plot
pdf("Q3.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of error ratio", 
     main = expression(Q[3]), xaxt = "n", cex.main = 1.4, cex.lab = 1.4, yaxt = "n")
lines(q, saddle_relative, cex = 1.2, lty = 4)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.2 %*% 10^-1), expression(1.2 %*% 10^-3), 
                                                  expression(1 %*% 10^-5), expression(1.3 %*% 10^-7), expression(1.3 %*% 
                                                                                                                                                                   10^-9), expression(1.3 %*% 10^-11), expression(1.4 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite",
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()


## Q4
## read non-zero eigenvalues
eigen = readRDS("eigenVal7000.rds")
q = (c(seq(8, by = 1.5, length.out = 4), 14.5, 18, 21.3, 24.5, 28) * 10^5)

## Farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 3985), eps = 1e-16, maxit = 2.14e+09)$Qq))

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 3985), acc = 1e-16, lim = 5e+07)$Qq))

## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 3985))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 3985), a = eigen, method = "satterthwaite"))

# saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 3985), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## plot
pdf("Rplots1.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of error ratio", 
     main = expression(Q[4]), xaxt = "n", cex.main = 1.4, cex.lab = 1.4, yaxt = "n")
lines(q, saddle_relative, cex = 1.2, lty = 4)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.8 %*% 10^-1), expression(1.9 %*% 10^-3), 
                                                  expression(1.3 %*% 10^-5), expression(1.1 %*% 10^-7), expression(1.2 %*% 
                                                                                                                                                                     10^-9), expression(1.6 %*% 10^-11), expression(1.4 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite", 
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()



## Q5
## read non-zero eigenvalues
eigen = readRDS("eigenVal9000.rds")
q = c(12, 14, 16, 18, 20.5, 25, 29.5, 34.5, 39) * 10^5

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 4984), acc = 1e-16, lim = 5e+07)$Qq))

## farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 4984), eps = 1e-16, maxit = 2.14e+09)$Qq))

## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 4984))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 4984), a = eigen, method = "satterthwaite"))

## Saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 4984), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## plot
pdf("Q5.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of relative error", 
     main = expression(Q[5]), xaxt = "n", yaxt = "n", cex.main = 1.4, cex.lab = 1.4)
lines(q, saddle_relative, lty = 4, cex = 1.2)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.1 %*% 10^-1), expression(1.4 %*% 10^-3), 
                                                  expression(1.4 %*% 10^-5), expression(1.5 %*% 10^-7), expression(1.8 %*% 
                                                                                                                                                                     10^-9), expression(1.5 %*% 10^-11), expression(1.8 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite", 
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()


## Q6
## read non-zero eigenvalues
eigen = readRDS("eigenVal20000.rds")
q = c(seq(6.1, by = 1, length.out = 3), 9.3, 10.5, 12.9, 15.4, 17.9, 20.3) * 10^6

## Davies
(dav = sapply(q, function(x) davies(x, eigen, rep(1, 11259), acc = 1e-16, lim = 5e+07)$Qq))
## Farebrother
(far = sapply(q, function(x) farebrother(x, eigen, rep(1, 11259), eps = 1e-16, maxit = 2.14e+09)$Qq))
## Liu
(l = sapply(q, function(x) liu(x, eigen, rep(1, 11259))))

## Satterthwaite
(satter = pchisqsum(q, lower.tail = FALSE, df = rep(1, 11259), a = eigen, method = "satterthwaite"))

## saddlepoint
(saddle = pchisqsum(q, lower.tail = FALSE, df = rep(1, 11259), a = eigen, method = "sad"))

## log of error ratio
l_relative = log10(l/dav)
satter_relative = log10(satter/dav)
saddle_relative = log10(saddle/dav)

## plot
pdf("Q6.pdf", height = 5, width = 7)
plot(q, saddle_relative, ylim = c(-2, 0.15), xlab = "Corresponding p-value", ylab = "Log of error ratio", 
     main = expression(Q[6]), xaxt = "n", yaxt = "n", cex.main = 1.4, cex.lab = 1.4)
lines(q, saddle_relative, lty = 4, cex = 1.2)
points(q, satter_relative, col = 3, cex = 1.2, pch = 2)
lines(q, satter_relative, col = 3, cex = 1.2, lty = 4)
points(q, l_relative, col = 2, cex = 1.2, pch = 5)
lines(q, l_relative, col = 2, cex = 1.2, lty = 4)
axis(1, at = q[c(1, 3, 5, 6, 7, 8, 9)], label = c(expression(1.5 %*% 10^-1), expression(1.8 %*% 10^-3), 
                                                  expression(1.6 %*% 10^-5), expression(1.6 %*% 10^-7), expression(1.3 %*% 
                                                                                                                                                                     10^-9), expression(1.2 %*% 10^-11), expression(1.2 %*% 10^-13)), cex.axis = 1)
axis(2, at = c(-2, -1.5, -1, -0.5, 0), labels = c("-2.0", "-1.5", "-1.0", "-0.5", 
                                                  "0.0"), las = 3, cex.axis = 1.2)
legend("bottomright", col = c(1, 3, 2), cex = 1.2, pch = c(1, 2, 5), lty = 4, legend = c("Saddlepoint", "Satterthwaite", 
                                                                                         "Liu-Tang-Zhang"), bty = "n")
dev.off()