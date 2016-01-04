### R code from vignette source 'CompareGroups.Rnw'

###################################################
### code chunk number 1: CompareGroups.Rnw:29-35
###################################################
# Load the stats, smwrData, and smwrStats packages
library(stats)
library(smwrData)
library(smwrStats)
# Get the dataset
data(AppalachianSpecCap)


###################################################
### code chunk number 2: CompareGroups.Rnw:43-51
###################################################
# setSweave is a specialized function that sets up the graphics page for
# Sweave scripts. For interactive use, it should be removed and the
# default setting for set.up can be used.
setSweave("CG01", 5, 5)
with(AppalachianSpecCap, boxPlot(LogSpecCap, group=RockType,
																 Box=list(type="tukey")))
# Required call to close PDF output graphics
graphics.off()


###################################################
### code chunk number 3: CompareGroups.Rnw:53-55
###################################################
cat("\\includegraphics{CG01.pdf}\n")
cat("\\paragraph{}\n")


###################################################
### code chunk number 4: CompareGroups.Rnw:65-67
###################################################
# Perform the analysis using a formula
kruskal.test(LogSpecCap ~ RockType, data=AppalachianSpecCap)


###################################################
### code chunk number 5: CompareGroups.Rnw:72-76
###################################################
# Perform the MCT using the default Tukey method for determining the
# critical value for separating groups.
with(AppalachianSpecCap, multicomp.test(LogSpecCap, RockType,
																				method="non"))


###################################################
### code chunk number 6: CompareGroups.Rnw:86-88
###################################################
# Perform the analysis using a formula
oneway.test(LogSpecCap ~ RockType, data=AppalachianSpecCap)


###################################################
### code chunk number 7: CompareGroups.Rnw:93-97
###################################################
# Perform the MCT using the default Tukey method for determining the
# critical value for separating groups.
with(AppalachianSpecCap, multicomp.test(LogSpecCap, RockType,
																				method="para"))


