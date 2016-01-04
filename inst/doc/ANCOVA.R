### R code from vignette source 'ANCOVA.Rnw'

###################################################
### code chunk number 1: ANCOVA.Rnw:29-35
###################################################
# Load the smwrStats and smwrData packages
library(smwrStats)
library(smwrData)
# Get the dataset
data(UraniumTDS)
head(UraniumTDS)


###################################################
### code chunk number 2: ANCOVA.Rnw:45-50
###################################################
# Create the ANCOVA model
UTDS.anc <- lm(log(Uranium) ~ HCO3*log(TDS), data=UraniumTDS)
# Perform the diagnostics after seleting the "best" subset
# The trace can be instructive, but is not necessary.
UTDS.best <- ancovaReg(UTDS.anc, trace=TRUE)


###################################################
### code chunk number 3: ANCOVA.Rnw:60-62
###################################################
# Print the final ANCOVA model
print(UTDS.best)


###################################################
### code chunk number 4: ANCOVA.Rnw:74-82
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("ancplot01", 5, 5)
plot(UTDS.best, which=1, set.up=FALSE)
# Required call to close PDF output graphics

graphics.off()


###################################################
### code chunk number 5: ANCOVA.Rnw:84-86
###################################################
cat("\\includegraphics{ancplot01.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 6: ANCOVA.Rnw:91-99
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("ancplot02", 5, 5)
plot(UTDS.best, which=2, set.up=FALSE)
# Required call to close PDF output graphics

graphics.off()


###################################################
### code chunk number 7: ANCOVA.Rnw:101-103
###################################################
cat("\\includegraphics{ancplot02.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 8: ANCOVA.Rnw:114-121
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("ancplot03", 5, 5)
plot(UTDS.best, which=3, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 9: ANCOVA.Rnw:123-125
###################################################
cat("\\includegraphics{ancplot03.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 10: ANCOVA.Rnw:135-142
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("ancplot04", 5, 5)
plot(UTDS.best, which=5, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 11: ANCOVA.Rnw:144-146
###################################################
cat("\\includegraphics{ancplot04.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 12: ANCOVA.Rnw:156-163
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("ancplot05", 5, 5)
plot(UTDS.best, which=6, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 13: ANCOVA.Rnw:165-167
###################################################
cat("\\includegraphics{ancplot05.pdf}\n")
cat("\\paragraph{}\n")


