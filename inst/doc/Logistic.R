### R code from vignette source 'Logistic.Rnw'

###################################################
### code chunk number 1: Logistic.Rnw:23-29
###################################################
# Load the smwrStats and smwrData packages
library(smwrStats)
library(smwrData)
# Get the dataset
data(PugetNitrate)
head(PugetNitrate)


###################################################
### code chunk number 2: Logistic.Rnw:40-45
###################################################
# Create the logistic regression model
PSNO3.1 <- glm(formula = nitrate >= 3 ~ wellmet, family = binomial, 
    data = PugetNitrate,  na.action = na.exclude)
# Print the summary 
print(summary(PSNO3.1))


###################################################
### code chunk number 3: Logistic.Rnw:52-73
###################################################
# Run the H-L test
PSNO3.1.hl <- hosmerLemeshow.test(PSNO3.1)
print(PSNO3.1.hl)
# Added fitted values to dataset for line in figure 2, and order
PugetNitrate$fits <- fitted(PSNO3.1)
OrderFits <- order(PugetNitrate$fits)
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("binplot01", 5, 5)
with(PugetNitrate, xyPlot(wellmet[OrderFits], fits[OrderFits],
    Plot=list(what="lines"),
    xaxis.range=c(0, 200),
    yaxis.range=c(0, .25),
    xtitle="Well Depth, in meters",
    ytitle="Estimated Pobability"))
# Add the observed frequencies
with(PSNO3.1.hl$estimate, addXY(wellmet, Counts/Size,
    Plot=list(what="points")))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 4: Logistic.Rnw:75-77
###################################################
cat("\\includegraphics{binplot01.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 5: Logistic.Rnw:84-86
###################################################
# Run the H-L test with 12 groups
hosmerLemeshow.test(PSNO3.1, 12)


###################################################
### code chunk number 6: Logistic.Rnw:91-93
###################################################
# Compute the area under the ROC
roc(PSNO3.1)


###################################################
### code chunk number 7: Logistic.Rnw:102-109
###################################################
# Create the logistic regression model
PSNO3.3 <- glm(formula = nitrate >= 3 ~ wellmet + l20 + l10, 
    family = binomial, subset = surfgeo == "Coarse",
    data = PugetNitrate,  na.action = na.omit)
# Create the assessment and print it
PSNO3.3.br <- binaryReg(PSNO3.3)
print(PSNO3.3.br)


