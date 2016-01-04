library(smwrStats)
# Set parameters for simulation
Seed <- 12319L
# Enter the correlation between X and Y
Cor <- 0.8
# Enter the sample size
Ns <- 25
# Enter the number of MC repetitions
NB <- 100
# Enter the critcal p-value for correlation, one-sided
Pc <- 0.05
# Enter the prediction probability level
PP <- 0.1
# Enter the ratio of the variances of y and x
Vrat <- 1.5
# Enter the offset to Y (added to Y, affects log-normal values)
Yoff <- 1.0
#-----------------------------------
## Create the population of Npop cases, max(2000, 100/PP)
Npop <- as.integer(max(2000, 100/PP))
set.seed(Seed)
Mat <- matrix(rnorm(Npop * 2L), ncol=2) %*% 
  chol(matrix(c(1, Cor*sqrt(Vrat), Cor*sqrt(Vrat), Vrat), ncol=2))
## Create storage vectors for the results
mvB <- double(NB) # The coefficient
crV <- double(NB) # The correlation
crP <- double(NB) # The attained p-value of the corr
prV <- double(NB) # The predicted Values
jkV <- double(NB) # The jackknife variance
mdV <- double(NB) # The move.1 variance 
jklV <- double(NB) # The jackknife variance of log normal
mdlV <- double(NB) # The move.1 variance of log normal
## Known, true values
prX <- quantile(Mat[, 1L], probs=PP) # Actual population values!
prY <- quantile(Mat[, 2L], probs=PP) + Yoff
trB <- sd(Mat[, 2L])/sd(Mat[, 1L])
COR <- cor(Mat[, 1L], Mat[, 2L])
## Do the MC simulations
for(i in seq(NB)) {
 rsamp <- sample(Npop, Ns)
 X <- Mat[rsamp, 1L]
 Y <- Mat[rsamp, 2L] + Yoff
 Mv <- move.1(Y ~ X)
 mvB[i] <- Mv$coefficients[2L]
 crV[i] <- Mv$R
 crP[i] <- 1 - pt(sqrt(Ns-2) * crV[i]/(1 - crV[i]^2), Ns - 2L)
 # normal 
 Mv.jk <- jackknifeMove.1(Mv, newdata=data.frame(X=prX))
 prV[i] <- Mv.jk[1L, 1L] # same as from predict
 jkV[i] <- Mv.jk[1L, 3L]
 mdV[i] <- predict(Mv, data.frame(X=prX), var.fit=TRUE)[,2L]
 # log normal
 X2 <- exp(X)
 Y2 <- exp(Y)
 MvX <- move.1(Y2 ~ X2, distribution="log")
 jklV[i] <- jackknifeMove.1(MvX, newdata=data.frame(X2=exp(prX)))[3L]
 mdlV[i] <- predict(MvX, data.frame(X2=exp(prX)), var.fit=TRUE)[,2L]
}
## The results:
## How reliable is the sign of the correlation?
Nexc <- sum(crP >= Pc)
Nnegr <- sum(crV < 0)
cat("Number of excluded cases based on critical P: ", Nexc,
    ", number of negative R: ", Nnegr, "\n", sep="")

## Compute the statistics for the coefficient:
MnB <- mean(mvB)
MnCB <- mean(mvB[crP < Pc])
cat("\nTrue Beta: ", round(trB, 4), ",\n Mean B: ", round(MnB, 4),
    ",\n Conditional Mean B: ", round(MnCB, 4), "\n", sep="")

VarB <- trB^2*(1 - COR^2)/(Ns - 2)
mseB <- mean((mvB - trB)^2)
mseCB <- mean((mvB[crP < Pc] - trB)^2)
cat("\nVariance of Beta: ", round(VarB, 4),
    ",\n MSE B: ", round(mseB, 4),
    ",\n Conditional MSE B: ", round(mseCB, 4), "\n", sep="")

## Compute the statistics for the prediction:
MnP <- mean(prV)
MnCP <- mean(prV[crP < Pc])
cat("\nTrue predicted value: ", round(prY, 4), 
    ",\n mean pred.: ", round(MnP, 4),
    ",\n cond. mean pred.: ", round(MnCP, 4), "\n", sep="")

mseP <- mean((prV - prY)^2)
mseCP <- mean((prV[crP < Pc] - prY)^2)
cat("\nVariance pred.:", round(mean(mdV), 4),
		",\n Jackknife mean variance: ", round(mean(jkV), 4),
    ",\n Bootstrap MSE pred.: ", round(mseP, 4),
    ",\n conditional MSE pred.: ", round(mseCP, 4), "\n", sep="")

## And for the log-normal prediction
mselP <- mean((exp(prV) - exp(prY))^2)
mselCP <- mean((exp(prV[crP < Pc]) - exp(prY))^2)
cat("\nMean lognormal prediction: ", round(mean(exp(prV)), 4),
		"\n Mean jackknife variance: ", round(mean(jklV), 4),
		"\n Mean Log-normal Expansion: ", round(mean(mdlV), 4),
		",\n Bootstrap MSE pred.: ", round(mselP, 4),
		",\n conditional MSE pred.: ", round(mselCP, 4),
		"\n")
