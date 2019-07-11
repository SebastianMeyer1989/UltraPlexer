args = commandArgs(trailingOnly=TRUE)

fn_input <- "C:\\Users\\diltheyat\\Documents\\Oxford\\documents\\analysis\\17 Mai 2018\\SIM1.classification_k19.evaluationAgainstMapping"

if(length(args) > 0)
{
	fn_input <- args[[1]]
}
if(nchar(fn_input) == 0)
{
	cat("Please provide evaluation file as first argument, e.g. Rscript plotCorrectness.R SIM1.classification_k19.evaluationAgainstMapping")
	
}	
stopifnot(nchar(fn_input) > 0)

fn_output <- paste(fn_input, ".pdf", sep = "")
pdf(fn_output, height = 10, width = 15)
par(mfrow=c(2,3)) 
D <- read.delim(fn_input, header = T, stringsAsFactors = F)
D <- D[D[["readLength"]] < 50000,]
D[["rL_2000"]] <- round(D[["readLength"]] / 2000) * 2000
D[["roundedIdty"]] <- round(D[["maximumIdentity"]]*100)
for(m in unique(D[["method"]]))
{
	D_m <- D[D[["method"]] == m, ]
	x_values <- c()

	# confidence_interval_sizes <- c()
	
	
	list_for_boxplot_abs <- list()
	list_for_boxplot_rel <- list()

	vector_for_readLength <- c()
	
	vector_for_barplot_001 <- c()
	vector_for_barplot_0001 <- c()
	
	
	vector_for_barplot_001_confidences <- c()
	vector_for_barplot_0001_confidences <- c()
			
	sorted_length_bins <- sort(unique(D_m[["rL_2000"]]))
	for(xv in sorted_length_bins)
	{
		x_values <- c(x_values, xv)
		y_values_abs <- D_m[["distanceToTruthAbs"]][D_m[["rL_2000"]] == xv]
		y_values_rel <- D_m[["distanceToTruthRel"]][D_m[["rL_2000"]] == xv]
		
		list_for_boxplot_abs[[as.character(xv)]] <- y_values_abs
		list_for_boxplot_rel[[as.character(xv)]] <- y_values_rel
		
		n_001 <- sum(y_values_rel < 0.01)
		n_0001 <- sum(y_values_rel < 0.001)
		
		n <- length(y_values_rel)
		p_001 <- n_001/n
		p_0001 <- n_0001/n
		
		vector_for_readLength <- c(vector_for_readLength, n)
		
		p_001_confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p_001*(1-p_001))
		p_0001_confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p_0001*(1-p_0001))
		
		# y_values <- c(y_values, p)
		
		# cat(m, xv, length(D_m[["correct"]][D_m[["rL_2000"]] == xv]), "\n")

		# n <- length(D_m[["correct"]][D_m[["rL_2000"]] == xv])
		# p <- mean(D_m[["correct"]][D_m[["rL_2000"]] == xv])
		
		# confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p*(1-p))
		
		#cat(m, xv, length(D_m[["correct"]][D_m[["rL_2000"]] == xv]), "\n")
		# y_values <- c(y_values, p)
		#confidence_interval_sizes <- c(confidence_interval_sizes, confidence_interval_size)
		
		vector_for_barplot_001 <- c(vector_for_barplot_001, p_001)
		vector_for_barplot_0001 <- c(vector_for_barplot_0001, p_0001)
		vector_for_barplot_001_confidences <- c(vector_for_barplot_001_confidences, p_001_confidence_interval_size)
		vector_for_barplot_0001_confidences <- c(vector_for_barplot_0001_confidences, p_0001_confidence_interval_size)		
	}
	

	x_values_idty <- c()
	vector_idty_for_barplot_0001 <- c()
	vector_idty_for_barplot_0001_confidences <- c()

	sorted_idty <- sort(unique(D_m[["roundedIdty"]]))
	for(xv in sorted_idty)
	{
		x_values_idty <- c(x_values_idty, xv)
		y_values_abs <- D_m[["distanceToTruthAbs"]][D_m[["roundedIdty"]] == xv]
		y_values_rel <- D_m[["distanceToTruthRel"]][D_m[["roundedIdty"]] == xv]

		n_001 <- sum(y_values_rel < 0.01)
		n_0001 <- sum(y_values_rel < 0.001)
		
		n <- length(y_values_rel)
		p_001 <- n_001/n
		p_0001 <- n_0001/n
				
		p_001_confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p_001*(1-p_001))
		p_0001_confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p_0001*(1-p_0001))
		
		# y_values <- c(y_values, p)
		
		cat(m, xv, length(D_m[["correct"]][D_m[["rL_2000"]] == xv]), "\n")
		vector_idty_for_barplot_0001 <- c(vector_idty_for_barplot_0001, p_001)
		vector_idty_for_barplot_0001_confidences <- c(vector_idty_for_barplot_0001_confidences, p_001_confidence_interval_size)		
	}
	
#	plot(x_values, y_values, xlab = "Read length bin", ylab = "Proportion correctly assigned", main = m)
#	for(i in 1:length(y_values))
#	{
		#lines(rep(x_values[[i]], 2), y_values[[i]] + (c(-1, 1) * confidence_interval_sizes[[i]]))
#	}

	hist(D_m[["distanceToTruthRel"]], nclass = 10000, xlim = c(0, 0.01), main = paste("Relative difference (mismatching bases / read length) - showing ", sprintf("%.2f", 100 * sum(D_m[["distanceToTruthRel"]] <= 0.01)/length(D_m[["distanceToTruthRel"]])), "% of reads.", sep = ""), xlab = "Relative difference" )
	
	plot(sorted_length_bins, vector_for_readLength, main = "Read length counts", xlab = "Read length")
	
	#boxplot(list_for_boxplot_abs, main = paste("Absolute difference (mismatching bases)", sep = ""), cex.axis = 1, ylim = c(min(c(list_for_boxplot_abs, recursive = T)), max(c(list_for_boxplot_abs, recursive = T))))
	boxplot(list_for_boxplot_rel, main = paste("Relative difference (mismatching bases / read length)", sep = ""), cex.axis = 1, ylim = c(-1, 1), xlab = "Read length")
	
	names(list_for_boxplot_abs) <- sorted_length_bins
	names(list_for_boxplot_rel) <- list_for_boxplot_rel
	
	#x_coords_001 <- barplot(vector_for_barplot_001, main = paste("Proportion reads with rel. diff. <1%", sep = ""))
	plot(sorted_length_bins, vector_for_barplot_001, main = paste("Proportion reads with rel. diff. <1%", sep = ""), ylim = c(0,1), xlab = "Read length")
	for(i in 1:length(vector_for_barplot_001_confidences))
	{
		# lines(rep(x_coords_001[[i]], 2), vector_for_barplot_001[[i]] + (c(-1, 1) * vector_for_barplot_001_confidences[[i]]), lwd = 2)
		lines(rep(sorted_length_bins[[i]], 2), vector_for_barplot_001[[i]] + (c(-1, 1) * vector_for_barplot_001_confidences[[i]]), lwd = 2)
	}	
	#x_coords_0001 <- barplot(vector_for_barplot_0001, main = paste("Proportion reads with rel. diff. <0.1%", sep = ""))
	plot(sorted_length_bins, vector_for_barplot_0001, main = paste("Proportion reads with rel. diff. <0.1%", sep = ""), ylim = c(0,1), xlab = "Read length")
	for(i in 1:length(vector_for_barplot_0001_confidences))
	{
		# lines(rep(x_coords_0001[[i]], 2), vector_for_barplot_0001[[i]] + (c(-1, 1) * vector_for_barplot_0001_confidences[[i]]), lwd = 2)
		lines(rep(sorted_length_bins[[i]], 2), vector_for_barplot_0001[[i]] + (c(-1, 1) * vector_for_barplot_0001_confidences[[i]]), lwd = 2)
	}	
	
	plot(sorted_idty, vector_idty_for_barplot_0001, main = paste("Proportion reads with rel. diff. <0.01% by identity", sep = ""), ylim = c(0,1), xlab = "Estimated identity")
	for(i in 1:length(vector_idty_for_barplot_0001_confidences))
	{
		# lines(rep(x_coords_0001[[i]], 2), vector_for_barplot_0001[[i]] + (c(-1, 1) * vector_for_barplot_0001_confidences[[i]]), lwd = 2)
		lines(rep(sorted_idty[[i]], 2), vector_idty_for_barplot_0001[[i]] + (c(-1, 1) * vector_idty_for_barplot_0001_confidences[[i]]), lwd = 2)
	}	
	
	mtext(paste("Classifier method: ", m), outer = TRUE, cex = 1.5, line = -2)	
}
dev.off()


	

stop()

fn_output <- paste(fn_input, ".pdf", sep = "")
pdf(fn_output)
D <- read.delim(fn_input, header = T, stringsAsFactors = F)
D[["rL_2000"]] <- round(D[["readLength"]] / 2000) * 2000
for(m in unique(D[["method"]]))
{
	D_m <- D[D[["method"]] == m, ]
	x_values <- c()
	y_values <- c()
	confidence_interval_sizes <- c()
	for(xv in sort(unique(D_m[["rL_2000"]])))
	{
		x_values <- c(x_values, xv)
		
		n <- length(D_m[["correct"]][D_m[["rL_2000"]] == xv])
		p <- mean(D_m[["correct"]][D_m[["rL_2000"]] == xv])
		
		confidence_interval_size <- qnorm(0.975)*sqrt((1/n)*p*(1-p))
		
		cat(m, xv, length(D_m[["correct"]][D_m[["rL_2000"]] == xv]), "\n")
		y_values <- c(y_values, p)
		confidence_interval_sizes <- c(confidence_interval_sizes, confidence_interval_size)
	}
	plot(x_values, y_values, xlab = "Read length bin", ylab = "Proportion correctly assigned", main = m)
	for(i in 1:length(y_values))
	{
		lines(rep(x_values[[i]], 2), y_values[[i]] + (c(-1, 1) * confidence_interval_sizes[[i]]))
	}
}
dev.off()

cat("Generated file: ", fn_output, "\n")
