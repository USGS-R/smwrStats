### R code from vignette source 'RecordExtension.Rnw'

###################################################
### code chunk number 1: RecordExtension.Rnw:29-39
###################################################
# Load the smwrStats and dataRetrieval packages
library(smwrStats)
library(dataRetrieval)
# Get the datasets and rename columns 
NFYB <- readNWISdv("05292704", parameterCd="00060", startDate="2000-10-01", 
  endDate="2001-09-30")
NFYB <- renameNWISColumns(NFYB)
YB <- readNWISdv("05293000", parameterCd="00060", startDate="2000-10-01", 
  endDate="2001-09-30")
YB <- renameNWISColumns(YB)


###################################################
### code chunk number 2: RecordExtension.Rnw:49-54
###################################################
# Merge the data
YBM <- merge(NFYB, YB, by="Date", all=TRUE, suffixes=c(".NFYB", ".YB"))
# Construct and print the model.
YBM.m2ln <- move.2(Flow.YB ~ Flow.NFYB, data=YBM, distribution="lognormal")
print(YBM.m2ln)


###################################################
### code chunk number 3: RecordExtension.Rnw:59-66
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("graph01", 6, 6)
plot(YBM.m2ln, which=2, set.up=FALSE)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 4: RecordExtension.Rnw:76-83
###################################################
# Predict all values
YBM$Pred.ln <- predict(YBM.m2ln)
setSweave("graph02", 6, 6)
AA.pl <- with(YBM, timePlot(Date, Pred.ln, yaxis.log=TRUE))
AA.pl <- with(YBM, addXY(Date, Flow.YB, 
  Plot=list(what="lines", color="green"), current=AA.pl))
graphics.off()


###################################################
### code chunk number 5: RecordExtension.Rnw:98-111
###################################################
# Construct and print the power transforms for multivariate normality
YBM.bc <- optimBoxCox(YBM[c("Flow.YB", "Flow.NFYB")])
print(YBM.bc)
# Construct and print the model.
YBM.m2bc <- move.2(Flow.YB ~ Flow.NFYB, data=YBM, distribution=YBM.bc)
print(YBM.m2bc)
# Predict all values
YBM$Pred.bc <- predict(YBM.m2bc)
setSweave("graph03", 6, 6)
AA.pl <- with(YBM, timePlot(Date, Pred.bc, yaxis.log=TRUE))
AA.pl <- with(YBM, addXY(Date, Flow.YB, 
  Plot=list(what="lines", color="green"), current=AA.pl))
graphics.off()


