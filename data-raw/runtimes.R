####
#
# Experimental runtimes for the five methods, for subsampled TCGA
#
####

# Greedy Forward Selection implementation to supplement method comparisons
source("greedy.R")

library(oscar) # OSCAR
library(glmnet) # LASSO
library(APML0) # APML0
library(ncvreg) # SCAD

# Example code for validating three transcriptomics data across multiple methods
library(curatedPCaData) # Version 0.9.42
library(survival) # Surv, coxph, etc

###
#
# The Cancer Genome Atlas (TCGA) PRAD data subset; used both for training (3/4) and held-out validation (1/4)
#
##
X1_samps <- rownames(colData(mae_tcga))[which(
	# Subset to primary samples
	colData(mae_tcga)$sample_type == "primary" & 
	# No NA recurrence statuses
	!is.na(colData(mae_tcga)$disease_specific_recurrence_status) & 
	# No NA recurrence follow-up times
	!is.na(colData(mae_tcga)$days_to_disease_specific_recurrence) &
	# Has to be present in transcriptomics
	rownames(colData(mae_tcga)) %in% colnames(mae_tcga[["gex.rsem.log"]])
	)]
Y1 <- Surv(event=colData(mae_tcga)[X1_samps, "disease_specific_recurrence_status"], time=colData(mae_tcga)[X1_samps, "days_to_disease_specific_recurrence"]+1)
X1 <- mae_tcga[["gex.rsem.log"]][,X1_samps]
X1 <- X1 |>
	(\(x) { x[apply(X1, MARGIN=1, FUN=\(q){ !any(table(q) > 0.5*length(q) ) }),] })() |> # Only include non-redundant variables; no more half of the samples can have the exactly same value
	(\(x) { x2 <- apply(x, MARGIN=1, FUN=scale); dimnames(x2) <- dimnames(t(x)); x2})() # z-score transformation to make datasets more comparable

# Create a grid of dimensionality (p) and sample sizes (n) to sample
runtimegrid <- expand.grid(
	# Dimensions
	p = c(50, 100, 500, 1000, 5000, 10000), 
	# Sample sizes
	n = c(50, 100, 200, 300), 
	# Methods
	method = c("oscar", "glmnet", "apml0", "ncvreg", "greedy")
)

# Store runtimes and run parameters
runs <- list()

for(i in 1:nrow(runtimegrid)){
	p <- runtimegrid$p[i]
	n <- runtimegrid$n[i]
	method <- runtimegrid$method[i]
	# Set a seet for comparability / reproducibility
	set.seed(p+n)
	cat(paste("method:", method, "\n"))
	cat(paste("p:", p, "\n"))
	cat(paste("n:", n, "\n"))
	cat("\n\n")
	# Down-sampled dimension (p) and sample size (n) from the whole TCGA transcriptomics data
	rows <- sample(1:nrow(X1), size=n)
	cols <- sample(1:ncol(X1), size=p)
	X <- X1[rows, cols]
	Y <- Y1[rows]
	if(method == "oscar"){
		time1 <- Sys.time() 
		fit <- oscar(x = X, y = Y, kmax = 50, family = "cox", in_selection=3, percentage=2*(10^-4), solver="LMBM");
		time2 <- Sys.time()
		runs[[length(runs)+1]] <- data.frame(time = time2 - time1, n = n, p = p, method = "oscar")
	}else if(method == "glmnet"){
		time1 <- Sys.time()
		fit <- glmnet(x = X, y = Y, family = "cox");
		time2 <- Sys.time()
		runs[[length(runs)+1]] <- data.frame(time = time2 - time1, n = n, p = p, method = "glmnet")				
	}else if(method == "apml0"){
		time1 <- Sys.time()
		APML0(x=X, y=Y, family="cox", iL0=TRUE, icutB = TRUE, keep.beta = T, nfolds = 1, foldid = NULL)
		time2 <- Sys.time()
		runs[[length(runs)+1]] <- data.frame(time = time2 - time1, n = n, p = p, method = "apml0")				
	}else if(method == "ncvreg"){
		time1 <- Sys.time()
		fit <- ncvsurv(X = X, y = Y, penalty = "SCAD")
		time2 <- Sys.time()
		runs[[length(runs)+1]] <- data.frame(time = time2 - time1, n = n, p = p, method = "ncvreg")				
	}else if(method == "greedy"){
		time1 <- Sys.time()
		fit <- greedyfw(x = X, y = Y, maxk = 50, seed = p+n, runCV = FALSE, verb = 0)
		time2 <- Sys.time()
		runs[[length(runs)+1]] <- data.frame(time = time2 - time1, n = n, p = p, method = "greedy")				
	}
	cat(paste("At length", length(runs), "of", nrow(runtimegrid), "\n"))
}

# Aggregate run times into a single matrix
runtimes <- do.call("rbind", runs)

runtimes$time <- as.numeric(runtimes$time)

library(ggplot2)

# Generate suitably paneled run time benchmark across methods
pdf("Fig_Runtimes.pdf", width=12, height=3.5)
ggplot(data = runtimes, aes(x = p, y = time, color = as.factor(n))) +
	geom_line(size = 1) +
	# Put time ticks at 0.1 seconds, 1 second, 5 seconds, 10 seconds, 60 seconds (1min), 600 seconds (10 minutes), 1800 seconds (30 minutes), 3600 seconds (1 hour)
	scale_y_continuous(trans='log10', breaks=c(0.1, 1, 5, 10, 60, 600, 1800, 3600)) +
	scale_color_manual(breaks = c("50", "100", "200", "300"), values = c("#7edaffff", "#00b0f6ff", "#51a7eaff", "#000000ff")) +
	facet_grid(. ~ factor(method, levels=c("oscar", "glmnet", "apml0", "ncvreg", "greedy")))
dev.off()
	
save.image("runtimes_workspace.RData")


	
	

