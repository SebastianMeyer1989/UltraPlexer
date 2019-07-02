fn_in <- "/gpfs/project/dilthey/projects/multiplexer/samples.txt.readStats"
D <- read.delim(fn_in, header = T)

table(D[["source"]])/length(D[["source"]])

sources <- unique(D[["source"]])
sources_c <- rainbow(length(sources))


D[["matchProp"]] <- -1
D[["rL_bin"]] <- round(D[["readLength"]]/2000)*2000

colNames_sources <- c()
for(S in sources)
{
	idx_S <- (D[["source"]] == S)
	D[["matchProp"]][idx_S] <- D[[paste("X", S, sep = "")]][idx_S]/D[["readLength"]][idx_S]
	colNames_sources <- c(colNames_sources, paste("X", S, sep = ""))
}
D[["bestNMatches"]] <- apply(D[,colNames_sources], 1, max)
D[["compatibleSamples"]]  <- 0
for(S in sources)
{
	D[["compatibleSamples"]] <- D[["compatibleSamples"]] + (D[[paste("X", S, sep = "")]] == D[["bestNMatches"]])
}

pdf(fn_out)

fn_out <- paste(fn_in, ".pdf", sep = "")

pdf(fn_out)
densities <- list()
densities_x <- c()
densities_y <- c()
for(sI in 1:length(sources))
{
	S <- sources[[sI]]
	S_d <- density(D[["readLength"]][D[["source"]] == S])
	densities[[sI]] <- S_d
	densities_x <- c(densities_x, S_d$x)
	densities_y <- c(densities_y, S_d$y)
}
plot(0, 0, col = "white", xlim = c(min(densities_x), max(densities_x)), ylim = c(min(densities_y), max(densities_y)))
for(sI in 1:length(sources))
{
	lines(densities[[sI]][["x"]], densities[[sI]][["y"]], col = sources_c[[sI]])
}
legend("topright", legend = sources, fill = sources_c)

matchprop <- list()
matchprop_xv_all <- c()
matchprop_yv_all <- c()
for(bin in sort(unique(D[["rL_bin"]])))
{
	idx2 <- ((D[["rL_bin"]] == bin))
	mean_matchprop <- mean(D[["matchProp"]][idx2])
	matchprop_xv_all <- c(matchprop_xv_all, bin)
	matchprop_yv_all <- c(matchprop_yv_all, mean_matchprop)
}
for(sI in 1:length(sources))
{
	S <- sources[[sI]]

	matchprop[[sI]] <- list()
	matchprop_xv <- c()
	matchprop_yv <- c()
	for(bin in sort(unique(D[["rL_bin"]])))
	{
		idx2 <- ((D[["rL_bin"]] == bin) & (D[["source"]] == S))
		mean_matchprop <- mean(D[["matchProp"]][idx2])
		matchprop_xv <- c(matchprop_xv, bin)
		matchprop_yv <- c(matchprop_yv, mean_matchprop)
	}
	matchprop[[sI]]$x <- matchprop_xv
	matchprop[[sI]]$y <- matchprop_yv
}
plot(0, 0, col = "white", xlim = c(0, max(D[["readLength"]])), ylim = c(0, 1))
for(sI in 1:length(sources))
{
	lines(matchprop[[sI]][["x"]], matchprop[[sI]][["y"]], col = sources_c[[sI]])
}
lines(matchprop_xv_all, matchprop_yv_all, col = "black", lwd = 2)
legend("topright", legend = sources, fill = sources_c)

compatibleSamples_xv <- c()
compatibleSamples_yv <- c()
for(bin in sort(unique(D[["rL_bin"]])))
{
	idx2 <- ((D[["rL_bin"]] == bin))
	mean_compatibleSamples <- mean(D[["compatibleSamples"]][idx2])
	compatibleSamples_xv <- c(compatibleSamples_xv, bin)
	compatibleSamples_yv <- c(compatibleSamples_yv, mean_compatibleSamples)
}
plot(compatibleSamples_xv, compatibleSamples_yv, main = "Mean compatible samples by read length", xlab = "Read length", ylab = "Mean compatible samples")

dev.off()

