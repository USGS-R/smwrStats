# Perform a simple seasonal wave analysis
library(smwrStats)
# Create the data set. These are Atrazine concentrations from the 
#  White River at Hazleton, IN. No data were censored in this time frame.
WHRatra <- data.frame(
CDATE=as.Date(c("1992-10-19", "1992-11-18", "1992-12-21", "1993-01-12",
"1993-02-18", "1993-03-30", "1993-04-12", "1993-04-26", "1993-05-03",
"1993-05-10", "1993-05-17", "1993-05-24", "1993-06-01", "1993-06-07",
"1993-06-14", "1993-06-21", "1993-06-28", "1993-07-06", "1993-07-12",
"1993-07-19", "1993-07-26", "1993-08-02", "1993-08-09", "1993-08-16",
"1993-08-23", "1993-08-30", "1993-09-07", "1993-09-20", "1993-10-12",
"1993-11-18", "1993-11-23", "1993-11-30", "1993-12-15", "1994-01-12",
"1994-02-15", "1994-03-05", "1994-04-11", "1994-05-09", "1994-05-23",
"1994-06-06", "1994-06-20", "1994-07-05", "1994-07-18", "1994-08-27",
"1994-09-12")),
CFLOW=c(4290, 31400, 7770, 39700, 9310, 20900, 31300, 18700, 16800, 18200,
15600, 9150, 6330, 10400, 9760, 12600, 9160, 21200, 11100, 11200, 6090, 4530,
3560, 6910, 23800, 6700, 22700, 7620, 13200, 87900, 113000, 30300, 19400,
11500, 13600, 13400, 25400, 33100, 11700, 5580, 4380, 6320, 4060, 1930, 1870),
Atra=c(0.280, 0.340, 0.180, 0.230, 0.100, 0.150, 0.540, 0.200, 1.100, 1.600,
1.300, 4.700, 1.900, 9.400, 12.800, 10.500, 7.800, 3.400, 2.900, 2.100, 1.100,
0.940, 0.670, 0.590, 0.450, 0.880, 0.260, 0.310, 0.210, 0.098, 0.190, 0.170,
0.140, 0.095, 0.060, 0.065, 0.740, 5.000, 4.000, 1.500, 3.900, 6.600, 1.800,
0.410, 0.430))

# Create a decimal time column
WHRatra$Dectime <- dectime(WHRatra$CDATE)
# The model is log(Atra) = B0 + B1*log(CFLOW) + B2*seasonalWave(Dectime) + E
# Fit the residuals from the simple regression between log(Atra) and log(CFLOW)
WHRatra$Resid <- residuals(lm(log(Atra) ~ log(CFLOW), data=WHRatra))
# First pass at the seasonal wave
WHR.sw <- with(WHRatra, seasonalPeak(Dectime, Resid))
# Confirm the seasonal wave:
# Select a single peak model--the second peak is too small to be modeled
#  and is best treated as noise.
# Click on the Red cross near the red vertical line in the graph and
#  enter 1 in the Console window.
WHR.swc <- confirm(WHR.sw)
print(WHR.swc)
# Now find the "best" model. In this case, we will not do an
#  exhaustive search because the timing of the peak was very distinct.
selBestWave(log(Atra) ~ log(CFLOW), data=WHRatra, dec.time="Dectime", wave.list=WHR.swc)
# The formula for model with the lowest AIC would be
# log(Atra) ~ log(CFLOW) + seasonalWave(Dectime, .4742, 3, 3)
