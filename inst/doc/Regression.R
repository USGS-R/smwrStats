### R code from vignette source 'Regression.Rnw'

###################################################
### code chunk number 1: Regression.Rnw:28-54
###################################################
# Load the smwrStats package
library(smwrStats)
# Create the Haan dataset
Haan1977 <- data.frame(
ROFF=c(17.38, 14.62, 15.48, 14.72, 18.37, 17.01, 18.2, 18.95, 13.94, 18.64,
 17.25, 17.48, 13.16),
PCIP=c(44.37, 44.09, 41.25, 45.5, 46.09, 49.12, 44.03, 48.71, 44.43, 47.72,
 48.38, 49, 47.03),
AREA=c(2.21, 2.53, 5.63, 1.55, 5.15, 2.14, 5.34, 7.47, 2.1, 3.89, 0.67,
 0.85, 1.72),
SLOPE=c(50, 7, 19, 6, 16, 26, 7, 11, 5, 18, 21, 23, 5),
LEN=c(2.38, 2.55, 3.11, 1.84, 4.14, 1.92, 4.73, 4.24, 2, 2.1, 1.15, 1.27,
 1.93),
PERIM=c(7.93, 7.65, 11.61, 5.31, 11.35, 5.89, 12.59, 12.33, 6.81, 9.87,
 3.93, 3.79, 5.19),
DI=c(0.91, 1.23, 2.11, 0.94, 1.63, 1.41, 1.3, 2.35, 1.19, 1.65, 0.62, 0.83,
 0.99),
Rs=c(0.38, 0.48, 0.57, 0.49, 0.39, 0.71, 0.27, 0.52, 0.53, 0.6, 0.48, 0.61,
 0.52),
FREQ=c(1.36, 2.37, 2.31, 3.87, 3.3, 1.87, 0.94, 1.2, 4.76, 3.08, 2.99,
 3.53, 2.33),
Rr=c(332, 55, 77, 68, 68, 230, 44, 72, 40, 115, 352, 300, 39)
)
# load the data library and get the Cuyahoga data
library(smwrData)
data(CuyahogaTDS)


###################################################
### code chunk number 2: Regression.Rnw:64-71
###################################################
# Create the allReg output dataset
HaanSub <- with(Haan1977, allReg(cbind(PCIP, AREA, SLOPE, LEN, PERIM, 
		DI, Rs, FREQ, Rr), ROFF, lin.dep=1))
# An alternative call, note the use of the drop argument
#HaanSub <- allReg(Haan1977[, -1], Haan1977[, 1, drop=FALSE], lin.dep=1)
# What are the "best" 5 models by Cp
head(HaanSub[order(HaanSub$Cp),])


###################################################
### code chunk number 3: Regression.Rnw:89-94
###################################################
# Create the regression model
Haan.lm <- lm(ROFF ~ PCIP + PERIM + Rr, data=Haan1977)
# Create the diagnostic object and print it.
Haan.reg <- multReg(Haan.lm)
print(Haan.reg)


###################################################
### code chunk number 4: Regression.Rnw:106-114
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("regplot01", 5, 5)
plot(Haan.reg, which=1, set.up=FALSE)
# Required call to close PDF output graphics

graphics.off()


###################################################
### code chunk number 5: Regression.Rnw:116-118
###################################################
cat("\\includegraphics{regplot01.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 6: Regression.Rnw:129-136
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("regplot02", 5, 5)
plot(Haan.reg, which=3, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 7: Regression.Rnw:138-140
###################################################
cat("\\includegraphics{regplot02.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 8: Regression.Rnw:150-157
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("regplot03", 5, 5)
plot(Haan.reg, which=5, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 9: Regression.Rnw:159-161
###################################################
cat("\\includegraphics{regplot03.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 10: Regression.Rnw:171-178
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("regplot04", 5, 5)
plot(Haan.reg, which=6, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 11: Regression.Rnw:180-182
###################################################
cat("\\includegraphics{regplot04.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 12: Regression.Rnw:192-199
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("regplot05", 5, 5)
plot(Haan.reg, which="PERIM", set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 13: Regression.Rnw:201-203
###################################################
cat("\\includegraphics{regplot05.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 14: Regression.Rnw:213-215
###################################################
# Compute the PRESS statistic
press(Haan.lm)


###################################################
### code chunk number 15: Regression.Rnw:220-224
###################################################
# The resdidual standard error is computed by the summary function.
summary(Haan.lm)
# But can easily be computed using rmse:
rmse(Haan.lm)


###################################################
### code chunk number 16: Regression.Rnw:229-231
###################################################
# The variance inflation factors:
vif(Haan.lm)


###################################################
### code chunk number 17: Regression.Rnw:239-244
###################################################
# Create the correlation structure, and print it:
Haan.cor <- cor.all(Haan1977)
print(Haan.cor, digits=3)
# Now summaryize the signficance of the realtions between ROFF and the other variables 
summary(Haan.cor, variable="ROFF")


###################################################
### code chunk number 18: Regression.Rnw:254-265
###################################################
# Create the regression model and print it:
TDS.lm <- lm(log(TDS) ~ log(Q) + fourier(TIME), data=CuyahogaTDS)
print(TDS.lm)
# The sum of the TDS data in the calibration dataset:
sum(CuyahogaTDS$TDS)
# The sum of the simple back-transformed predictions
sum(exp(predict(TDS.lm)))
# No the sume from each of the bais-corrected methods
sum(predictFerguson(TDS.lm))
sum(predictDuan(TDS.lm))
sum(predictMVUE(TDS.lm))


