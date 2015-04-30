### R code from vignette source 'Trends.Rnw'

###################################################
### code chunk number 1: Trends.Rnw:21-28
###################################################
# Load the stats, smwrData, and smwrStats packages
library(stats)
library(smwrData)
library(smwrStats)
# Get the datasets
data(ConecuhFlows)
data(KlamathTP)


###################################################
### code chunk number 2: Trends.Rnw:36-41
###################################################
setSweave("trend01", 5, 5)
with(ConecuhFlows, timePlot(Year, Flow, Plot=list(what="both")))
with(ConecuhFlows, refLine(horizontal=median(Flow)))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 3: Trends.Rnw:56-60
###################################################
# The default method, Wilcosxon rank-sum serial test
with(ConecuhFlows, serial.test(Flow))
# The runs test method
with(ConecuhFlows, serial.test(Flow, method="runs"))


###################################################
### code chunk number 4: Trends.Rnw:68-75
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("trend02", 5, 5)
with(KlamathTP, timePlot(sample_dt, TP_ss, Plot=list(what="points")))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 5: Trends.Rnw:84-88
###################################################
setSweave("trend03", 5, 5)
with(KlamathTP, xyPlot(Flow, TP_ss))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 6: Trends.Rnw:103-109
###################################################
# Construct the regular series: 96 monthly observations
KlamathTP.RS <- with(KlamathTP, regularSeries(TP_ss, sample_dt,
										 begin="1972-01-01", end="1980-01-01"))
# Do the analysis: Value is the name of the column in KlamathTP.RS
# that contains the TP data
with(KlamathTP.RS, seaken(Value, 12))


###################################################
### code chunk number 7: Trends.Rnw:116-127
###################################################
# Compute the flow-adjusted concentrations. Figure 2 serves as justification
# for the linear fit for these data.
KlamathTP$FAC <- residuals(lm(TP_ss ~ Flow, data=KlamathTP,
	na.action=na.exclude)) # required to preserve any missing values
# Construct the regular series: 96 monthly observations
KlamathTP.RS2 <- with(KlamathTP, regularSeries(FAC, sample_dt,
										 begin="1972-01-01", end="1980-01-01"))
# Do the analysis: Value is the name of the column in KlamathTP.RS
# that contains the FAC data
KTP <- with(KlamathTP.RS2, seaken(Value, 12))
print(KTP)


###################################################
### code chunk number 8: Trends.Rnw:142-150
###################################################
# Simple substitution for one left-censored value
KlamathTP <- transform(KlamathTP, TP_ss2 = ifelse(TP_rmk == "<", TP/2, TP),
											 Dectime = dectime(sample_dt))
# The trend analysis, residual plot review indicates that
# Flow must second-order 
KTP.lm <- lm(log(TP_ss2) ~ Dectime + quadratic(log(Flow)) + fourier(Dectime),
						 data=KlamathTP)
summary(KTP.lm)


