args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) > 1)
 
coverageFile <- args[[1]]
sampleID <- args[[2]]

histogramData <- read.delim(coverageFile, header = T, stringsAsFactors = F)

stopifnot((which.max(histogramData[[2]]) == 2) || (which.max(histogramData[[2]]) == 3) )
minimumValue <- histogramData[[2]][[2]]
minimumValueI <- 2
if(which.max(histogramData[[2]]) == 3)
{
	minimumValueI <- 3
}
while(histogramData[[2]][[minimumValueI+1]] <= minimumValue)
{
	minimumValueI <- minimumValueI + 1
	minimumValue <- histogramData[[2]][[minimumValueI]]
}

# cat(minimumValueI, "\n")
# cat(minimumValue, "\n")

nextPeakHeight <- max(histogramData[[2]][minimumValueI:length(histogramData[[2]])])
nextPeakHeight_x <- which.max(histogramData[[2]][minimumValueI:length(histogramData[[2]])])+minimumValueI

minimumValue_inBetween <- min(histogramData[[2]][minimumValueI:nextPeakHeight_x])
minimumValue_inBetween_x <- which.min(histogramData[[2]][minimumValueI:nextPeakHeight_x])+minimumValueI

if(minimumValue_inBetween_x > minimumValueI)
{
	minimumValueI <- minimumValue_inBetween_x
}

# cat("Resut: ", minimumValueI, "\n")

pdf(paste(coverageFile, ".pdf", sep = ""))
plot(histogramData[[1]], histogramData[[2]], main = paste("Coverage plot sample ", sampleID, "\n", "Threshold = ", minimumValueI, sep = ""), ylim = c(0, nextPeakHeight * 1.5), xlim = c(0, 2 * nextPeakHeight_x))
lines(c(minimumValueI, minimumValueI), c(0, nextPeakHeight), col = "red")
dev.off()

cat(minimumValueI, "\n", sep =  "", file = paste(coverageFile, "T", sep = ""))
