sampleFile <- "samples.txt"
classificationResults <- paste(sampleFile, '.classification', sep = "")
outputFile <- paste(classificationResults, '.plot.pdf', sep = "")
D <- read.delim(classificationResults, header = T)
truthValues <- unique(D[["truth"]])
stopifnot(length(truthValues) > 0)
pdf(outputFile)
for(tV in truthValues)
{
	tV_indices <- which(D[["truth"]] == tV)
	col_correct_p_1 <- paste("P", tV, "_proper" , sep = "")
	hist(D[[col_correct_p_1]][tV_indices], main = tV)
}
dev.off()