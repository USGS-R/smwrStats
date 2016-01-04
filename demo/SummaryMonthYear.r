# Example of computing summary statistics by month and year
library(lubridate) # needed for month and year functions
library(dataRetrieval) # Needed for data retrieval
# Get 2 water years of streamflow record for the Little Arkansas River
# near Sedgwick, Kansas
Lark <- renameNWISColumns(readNWISdv("07144100", "00060",
    startDate="2010-10-01", endDate="2012-09-30"))
# The interaction function can be used to combine factors or other grouping
# variables into a single grouping factor. It preserves the order and
# unused levels can be dropped if necessary.
Lark.sum <- aggregate(Flow ~ interaction(month(Date, label=TRUE),
    year(Date), drop=TRUE, sep=" "), data=Lark, FUN=mean)
# Change the name
names(Lark.sum)[1] <- "MonthYear"
# prnt the first few rows of the output
head(Lark.sum)
