### R code from vignette source 'LowFlowExtension.Rnw'

###################################################
### code chunk number 1: LowFlowExtension.Rnw:29-43
###################################################
# Load the smwrStats and dataRetrieval packages
library(smwrStats)
library(dataRetrieval)
# Get the reference datasets and Initial processing 
Passaic.Ref <- readNWISdv("01379000", parameterCd="00060",  startDate="1960-10-01", 
  endDate="2014-09-30")
Passaic.Ref <- renameNWISColumns(Passaic.Ref)
Whippany.Ref <- readNWISdv("01381500", parameterCd="00060",  startDate="1960-10-01",
  endDate="2014-09-30")
Whippany.Ref <- renameNWISColumns(Whippany.Ref)
# Get the partial record data
Passaic.PR <- readNWISmeas("01378700", endDate="2014-09-30")
# Need date only for matching
Passaic.PR <- transform(Passaic.PR, Date=as.Date(measurement_dt))


###################################################
### code chunk number 2: LowFlowExtension.Rnw:53-68
###################################################
# Compute the condition of the flow
Passaic.Ref <- transform(Passaic.Ref, Flow_cond=hysteresis(Flow))
Whippany.Ref <- transform(Whippany.Ref, Flow_cond=hysteresis(Flow))
# Determine base flow, will be TRUE for base flow
Passaic.Ref <- transform(Passaic.Ref, Base=Flow < median(Flow) &
  Flow_cond > -5 & Flow_cond <= 0)
# There are a few missing values in the flow record at 01381500
Whippany.Ref <- transform(Whippany.Ref, Base=Flow < median(Flow, na.rm=TRUE) & 
  Flow_cond > -5 & Flow_cond <= 0)
# Merge The reference flows with the measured flows
Passaic.Mrg <- merge(Passaic.PR, Passaic.Ref, by="Date")
Whippany.Mrg <- merge(Passaic.PR, Whippany.Ref, by="Date")
# Retain only the base-flow data
Passaic.Mrg <- subset(Passaic.Mrg, Base)
Whippany.Mrg <- subset(Whippany.Mrg, Base)


###################################################
### code chunk number 3: LowFlowExtension.Rnw:71-91
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("graph01", 6, 6)
AA.lo <- setLayout(num.rows=2)
# Plot the data
setGraph(1, AA.lo)
AA.pl <- with(Passaic.Mrg, xyPlot(Flow, discharge_va, xaxis.log=T, yaxis.log=T))
# Add the linear regresion line to asses the goodness of fit
addSLR(AA.pl)
addTitle("Passaic")
# Plot the data
setGraph(2, AA.lo)
AA.pl <- with(Whippany.Mrg, xyPlot(Flow, discharge_va, xaxis.log=T, yaxis.log=T))
addSLR(AA.pl)
addTitle("Whippany")
# Add the linear regresion line to asses the goodness of fit
addSLR(AA.pl)
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 4: LowFlowExtension.Rnw:106-113
###################################################
# Construct the and print the first model
Passaic1.mv1 <- move.1(discharge_va ~ Flow, data=Passaic.Mrg, distribution = "lognormal")
print(Passaic1.mv1)
# PLot the first diagnostic plot
setSweave("graph02", 6, 6)
plot(Passaic1.mv1, which=1, set.up=FALSE)
graphics.off()


###################################################
### code chunk number 5: LowFlowExtension.Rnw:123-130
###################################################
# Construct the and print the second model
Passaic2.mv1 <- move.1(discharge_va ~ Flow, data=Whippany.Mrg, distribution = "lognormal")
print(Passaic2.mv1)
# PLot the first diagnostic plot
setSweave("graph03", 6, 6)
plot(Passaic2.mv1, which=1, set.up=FALSE)
graphics.off()


###################################################
### code chunk number 6: LowFlowExtension.Rnw:140-144
###################################################
# Compute the 10th percentile of flow
Whippany.10 <- quantile(Whippany.Ref$Flow, probs=0.1, na.rm=TRUE)
# Use that to estiamte the 10th percentil at the partial record station
predict(Passaic2.mv1, newdata=data.frame(Flow=Whippany.10))


