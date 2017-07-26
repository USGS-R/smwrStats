### R code from vignette source 'Trends.Rnw'

###################################################
### code chunk number 1: Trends.Rnw:30-37
###################################################
# Load the stats, smwrData, and smwrStats packages
library(stats)
library(smwrData)
library(smwrStats)
# Get the datasets
data(ConecuhFlows)
data(KlamathTP)


###################################################
### code chunk number 2: Trends.Rnw:45-50
###################################################
setSweave("trend01", 5, 5)
with(ConecuhFlows, timePlot(Year, Flow, Plot=list(what="both")))
with(ConecuhFlows, refLine(horizontal=median(Flow)))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 3: Trends.Rnw:61-63
###################################################
# The Mann-Kendall trend test
with(ConecuhFlows, kensen.test(Flow, Year, 25))


###################################################
### code chunk number 4: Trends.Rnw:70-74
###################################################
# The default method, Wilcoxon rank-sum serial test
with(ConecuhFlows, serial.test(Flow))
# The runs test method
with(ConecuhFlows, serial.test(Flow, method="runs"))


###################################################
### code chunk number 5: Trends.Rnw:82-89
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("trend02", 5, 5)
with(KlamathTP, timePlot(sample_dt, TP_ss, Plot=list(what="points")))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 6: Trends.Rnw:98-102
###################################################
setSweave("trend03", 5, 5)
with(KlamathTP, xyPlot(Flow, TP_ss))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 7: Trends.Rnw:117-123
###################################################
# Construct the regular series: 96 monthly observations
KlamathTP.RS <- with(KlamathTP, regularSeries(TP_ss, sample_dt,
										 begin="1972-01-01", end="1980-01-01"))
# Do the analysis: Value is the name of the column in KlamathTP.RS
# that contains the TP data
with(KlamathTP.RS, seaken(Value, 12))


###################################################
### code chunk number 8: Trends.Rnw:130-141
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
### code chunk number 9: Trends.Rnw:156-164
###################################################
# Simple substitution for one left-censored value
KlamathTP <- transform(KlamathTP, TP_ss2 = ifelse(TP_rmk == "<", TP/2, TP),
											 Dectime = dectime(sample_dt))
# The trend analysis, residual plot review indicates that
# Flow must second-order 
KTP.lm <- lm(log(TP_ss2) ~ Dectime + quadratic(log(Flow)) + fourier(Dectime),
						 data=KlamathTP)
summary(KTP.lm)


###################################################
### code chunk number 10: Trends.Rnw:174-184
###################################################
# Get the data from the downloads folder for the report
RK3b <- read.delim("https://pubs.usgs.gov/sir/2005/5275/downloads/RK3b.txt",
  header=FALSE, skip=1)
names(RK3b) <- c("Year", "Site", "NH4")
# Reformat to a wide data frame group2row is in smwrBase and make the matrix
RK3b <- group2row(RK3b, "Year", "Site", "NH4")
RK3b.mat <- data.matrix(RK3b[, -1]) # column 1 is Year
# Preform the trend test
regken(RK3b.mat)
regken(RK3b.mat, correct=FALSE)


###################################################
### code chunk number 11: Trends.Rnw:189-191
###################################################
# Pretty print the Spearman correlation
printCor(cor(RK3b.mat, method="spearman", use="pair"), 0.5)


