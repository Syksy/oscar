# Greedy Forward Selection implementation to supplement method comparisons
source("greedy.R")

# Example code for validating three transcriptomics data across multiple methods
library(curatedPCaData) # Processing of PCa datasets
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
##
#
# Taylor et al. / MSKCC with recurrence and follow-up; use as external validation dataset
#
##
X2_samps <- rownames(colData(mae_taylor))[which(
	# Subset to primary samples
	colData(mae_taylor)$sample_type == "primary" & 
	# No NA recurrence statuses
	!is.na(colData(mae_taylor)$disease_specific_recurrence_status) & 
	# No NA recurrence follow-up times
	!is.na(colData(mae_taylor)$days_to_disease_specific_recurrence) &
	# Has to be present in transcriptomics
	rownames(colData(mae_taylor)) %in% colnames(mae_taylor[["gex.rma"]])
	)]
Y2 <- Surv(event=colData(mae_taylor)[X2_samps, "disease_specific_recurrence_status"], time=colData(mae_taylor)[X2_samps, "days_to_disease_specific_recurrence"])
X2 <- mae_taylor[["gex.rma"]][,X2_samps]
X2 <- X2 |>
	(\(x) { x2 <- apply(x, MARGIN=1, FUN=scale); dimnames(x2) <- dimnames(t(x)); x2})()  # z-score transformation to make datasets more comparable

##
#
# Sun et al. data matrix and binary recurrence outcome (no follow-up time available); use as external validation dataset
#
##
X3_samps <- rownames(colData(mae_sun))[which(
	# Subset to primary samples
	colData(mae_sun)$sample_type == "primary" & 
	# No NA recurrence statuses
	!is.na(colData(mae_sun)$disease_specific_recurrence_status) & 
	## No follow-up information! Run using logistic regression
	# No NA recurrence follow-up times
	#!is.na(colData(mae_sun)$days_to_disease_specific_recurrence) &
	# Has to be present in transcriptomics
	rownames(colData(mae_sun)) %in% colnames(mae_sun[["gex.rma"]])
	)]
Y3 <- as.numeric(colData(mae_sun)[X3_samps, "disease_specific_recurrence_status"])
X3 <- mae_sun[["gex.rma"]][,X3_samps]
X3 <- X3 |>
	(\(x) { x[apply(X1, MARGIN=1, FUN=\(q){ !any(table(q) > 0.5*length(q) ) }),] })() |> # Only include non-redundant variables; no more half of the samples can have the exactly same value
	(\(x) { x2 <- apply(x, MARGIN=1, FUN=scale); dimnames(x2) <- dimnames(t(x)); x2})()  # z-score transformation to make datasets more comparable

# Substitute difficult symbols with '.'
# e.g. coxph interface for greedy FS cannot handle variable name 'GENENAME-AS1" due to '-' being an operator
colnames(X1) <- gsub("-|/", ".", colnames(X1))
colnames(X2) <- gsub("-|/", ".", colnames(X2))
colnames(X3) <- gsub("-|/", ".", colnames(X3))

# Set aside one fourth of TCGA as internal validation
set.seed(123)
validation <- sample(1:nrow(X1), size=floor(nrow(X1)/4))
X4 <- X1[validation,]
Y4 <- Y1[validation]

X1 <- X1[-validation,]
Y1 <- Y1[-validation]

# Create an intersection of genes present in all studies (p=10,253)
commong <- intersect(intersect(colnames(X1), colnames(X2)), colnames(X3))

X1 <- X1[,commong] # TCGA training data
X2 <- X2[,commong] # Taylor et al.
X3 <- X3[,commong] # Sun et al.
X4 <- X4[,commong] # TCGA held-out data

#> dim(X1)
#[1]   303 10253
#> dim(X2)
#[1]   131 10253
#> dim(X3)
#[1]    79 10253
#> dim(X4)
#[1]   101 10253

# Required libraries for fits and validation
library(survival)
library(survminer)
library(pROC)

library(oscar)# OSCAR
library(glmnet) # glmnet
library(ncvreg) # SCAD
library(APML0) # APML0
# Greedy implementation defined above

## Fit OSCAR
set.seed(1)
fit_oscar <- oscar(x=X1, y=Y1, family="cox", kmax=50, in_selection=3, percentage=2*(10^-4), solver="LMBM")

## Fit Greedy FS
set.seed(1)
fit_greedy <- greedyfw(x = X1, y = Y1, verb = 1, maxk = 50, seed = 1, runCV = FALSE)

set.seed(1)
## glmnet (LASSO) fits
fit_glmnet <- glmnet(y = Y1, x = X1, family="cox")
cv.fit_glmnet <- cv.glmnet(y = Y1, x = X1, family="cox")

set.seed(1)
# APML0 fits
fit_apml0<- APML0(x=X1,y=Y1, family="cox", iL0=TRUE, icutB = TRUE, nfolds = 10, keep.beta = T)
bmat_apml0 <- as.matrix(fit_apml0$Beta)

set.seed(1)
# SCAD fits
fit_scad <- ncvsurv(X = X1, y = Y1, penalty = "SCAD")
cv.fit_scad <- cv.ncvsurv(X = X1, y = Y1, penalty = "SCAD")

# Transform Taylor et al. data matrix to comform with TCGA data matrix in case the variables are required
X2 <- X2 |> 
	(\(x){ apply(x, MARGIN=2, FUN=\(q){ q[is.na(q)] <- median(q, na.rm=TRUE); q }) })() |> # Median imputation if individual elements were missing from original data matrix
	(\(x){ x[,match(colnames(X1), colnames(x))]})() |> # Map to variables present in TCGA (X1)
	(\(x){ x[is.na(x)] <- median(x, na.rm=TRUE); x })() |> # Grand median imputed if NA columns were introduced when mapping to TCGA dimensions
	(\(x){ colnames(x) <- colnames(X1); x })() # Copy row names as originally presented in TCGA
	
# Transform Sun et al. data matrix to comform with TCGA data matrix in case the variables are required
X3 <- X3 |> 
	(\(x){ apply(x, MARGIN=2, FUN=\(q){ q[is.na(q)] <- median(q, na.rm=TRUE); q }) })() |> # Median imputation if individual elements were missing from original data matrix
	(\(x){ x[,match(colnames(X1), colnames(x))]})() |> # Map to variables present in TCGA (X1)
	(\(x){ x[is.na(x)] <- median(x, na.rm=TRUE); x })() |> # Grand median imputed if NA columns were introduced when mapping to TCGA dimensions
	(\(x){ colnames(x) <- colnames(X1); x })() # Copy row names as originally presented in TCGA

## OSCAR predictions	
	
# OSCAR predictions for Taylor et al.	
oscar_preds_X1 <- lapply(1:fit_oscar@kmax, FUN=\(k){
	predict(fit_oscar, newdata= X1, k=k, type="response")
})

# OSCAR predictions for Taylor et al.	
oscar_preds_X2 <- lapply(1:fit_oscar@kmax, FUN=\(k){
	predict(fit_oscar, newdata= X2, k=k, type="response")
})

# OSCAR predictions for Sun et al.	
oscar_preds_X3 <- lapply(1:fit_oscar@kmax, FUN=\(k){
	predict(fit_oscar, newdata= X3, k=k, type="response")
})

# OSCAR predictions for TCGA held out set	
oscar_preds_X4 <- lapply(1:fit_oscar@kmax, FUN=\(k){
	predict(fit_oscar, newdata= X4, k=k, type="response")
})


## Greedy predictions
	
# Greedy predictions for Taylor et al.	
greedy_preds_X1 <- lapply(1:50, FUN=\(k){
	# Coefs from TCGA
	model <- coxph(as.formula(paste("Y1 ~ ", paste(fit_greedy$selected[1:k], collapse=" +"))), data = as.data.frame(X1))
	# Predict to new data
	predict(model, newdata = as.data.frame(X1), type="lp")
})

# Greedy predictions for Taylor et al.	
greedy_preds_X2 <- lapply(1:50, FUN=\(k){
	# Coefs from TCGA
	model <- coxph(as.formula(paste("Y1 ~ ", paste(fit_greedy$selected[1:k], collapse=" +"))), data = as.data.frame(X1))
	# Predict to new data
	predict(model, newdata = as.data.frame(X2), type="lp")
})

# Greedy predictions for Sun et al.
greedy_preds_X3 <- lapply(1:50, FUN=\(k){
	# Coefs from TCGA
	model <- coxph(as.formula(paste("Y1 ~ ", paste(fit_greedy$selected[1:k], collapse=" +"))), data = as.data.frame(X1))
	# Predict to new data
	predict(model, newdata = as.data.frame(X3), type="lp")
})

# Greedy predictions for TCGA held out
greedy_preds_X4 <- lapply(1:50, FUN=\(k){
	# Coefs from TCGA
	model <- coxph(as.formula(paste("Y1 ~ ", paste(fit_greedy$selected[1:k], collapse=" +"))), data = as.data.frame(X1))
	# Predict to new data
	predict(model, newdata = as.data.frame(X4), type="lp")
})


## glmnet predictions

# glmnet predictions for Taylor et al.	
glmnet_preds_X1 <- lapply(1:ncol(fit_glmnet$beta), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- t(fit_glmnet$beta[,b])
	as.vector(beta %*% t(X1))
})

# glmnet predictions for Taylor et al.	
glmnet_preds_X2 <- lapply(1:ncol(fit_glmnet$beta), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- t(fit_glmnet$beta[,b])
	as.vector(beta %*% t(X2))
})

# glmnet predictions for Sun et al.	
glmnet_preds_X3 <- lapply(1:ncol(fit_glmnet$beta), FUN=\(b){
	# Multiply beta vector by X3 and extract the linear product per each lambda
	beta <- t(fit_glmnet$beta[,b])
	as.vector(beta %*% t(X3))
})

# glmnet predictions for TCGA held out
glmnet_preds_X4 <- lapply(1:ncol(fit_glmnet$beta), FUN=\(b){
	# Multiply beta vector by X4 and extract the linear product per each lambda
	beta <- t(fit_glmnet$beta[,b])
	as.vector(beta %*% t(X4))
})

## SCAD predictions

# SCAD predictions for Taylor et al.	
scad_preds_X1 <- lapply(1:ncol(fit_scad$beta), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- t(fit_scad$beta[,b])
	as.vector(beta %*% t(X1))
})

# SCAD predictions for Taylor et al.	
scad_preds_X2 <- lapply(1:ncol(fit_scad$beta), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- t(fit_scad$beta[,b])
	as.vector(beta %*% t(X2))
})

# SCAD predictions for Sun et al.	
scad_preds_X3 <- lapply(1:ncol(fit_scad$beta), FUN=\(b){
	# Multiply beta vector by X3 and extract the linear product per each lambda
	beta <- t(fit_scad$beta[,b])
	as.vector(beta %*% t(X3))
})

# SCAD predictions for TCGA held out
scad_preds_X4 <- lapply(1:ncol(fit_scad$beta), FUN=\(b){
	# Multiply beta vector by X4 and extract the linear product per each lambda
	beta <- t(fit_scad$beta[,b])
	as.vector(beta %*% t(X4))
})

## APML0 predictions

# APML0 predictions for Taylor et al.	
apml0_preds_X1 <- lapply(1:ncol(bmat_apml0), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- bmat_apml0[,b]
	as.vector(beta %*% t(X1))
})

# APML0 predictions for Taylor et al.	
apml0_preds_X2 <- lapply(1:ncol(bmat_apml0), FUN=\(b){
	# Multiply beta vector by X2 and extract the linear product per each lambda
	beta <- bmat_apml0[,b]
	as.vector(beta %*% t(X2))
})

# APML0 predictions for Sun et al.	
apml0_preds_X3 <- lapply(1:ncol(bmat_apml0), FUN=\(b){
	# Multiply beta vector by X3 and extract the linear product per each lambda
	beta <- bmat_apml0[,b]
	as.vector(beta %*% t(X3))
})

# APML0 predictions for TCGA held out
apml0_preds_X4 <- lapply(1:ncol(bmat_apml0), FUN=\(b){
	# Multiply beta vector by X3 and extract the linear product per each lambda
	beta <- bmat_apml0[,b]
	as.vector(beta %*% t(X4))
})



## OSCAR performance

# C-indices from OSCAR for predicting Taylor et al. for k=1:50
oscar_C_X1 <- unlist(lapply(oscar_preds_X1, FUN=\(x){ coxph(Y1 ~ x)[["concordance"]]["concordance"]  }))

# C-indices from OSCAR for predicting Taylor et al. for k=1:50
oscar_C_X2 <- unlist(lapply(oscar_preds_X2, FUN=\(x){ coxph(Y2 ~ x)[["concordance"]]["concordance"]  }))

# ROC-AUCs from OSCAR for predicting Sun et al. for k=1:50
oscar_ROCAUC_X3 <- unlist(lapply(oscar_preds_X3, FUN=\(x){ auc(roc(response=Y3, predictor=as.vector(x)))  }))

# C-indices from OSCAR for predicting Taylor et al. for k=1:50
oscar_C_X4 <- unlist(lapply(oscar_preds_X4, FUN=\(x){ coxph(Y4 ~ x)[["concordance"]]["concordance"]  }))


## Greedy FS performance

# C-indices
greedy_C_X1 <- unlist(lapply(greedy_preds_X1, FUN=\(x){ coxph(Y1 ~ x)[["concordance"]]["concordance"]  }))

# C-indices
greedy_C_X2 <- unlist(lapply(greedy_preds_X2, FUN=\(x){ coxph(Y2 ~ x)[["concordance"]]["concordance"]  }))

# ROC-AUCs
greedy_ROCAUC_X3 <- unlist(lapply(greedy_preds_X3, FUN=\(x){ auc(roc(response=Y3, predictor=x))  }))

# C-indices
greedy_C_X4 <- unlist(lapply(greedy_preds_X4, FUN=\(x){ coxph(Y4 ~ x)[["concordance"]]["concordance"]  }))


## glmnet performance

# C-indices
glmnet_C_X1 <- unlist(lapply(glmnet_preds_X1, FUN=\(x){ coxph(Y1 ~ x)[["concordance"]]["concordance"]  }))

# C-indices
glmnet_C_X2 <- unlist(lapply(glmnet_preds_X2, FUN=\(x){ coxph(Y2 ~ x)[["concordance"]]["concordance"]  }))

# ROC-AUCs
glmnet_ROCAUC_X3 <- unlist(lapply(glmnet_preds_X3, FUN=\(x){ auc(roc(response=Y3, predictor=x))  }))

# C-indices
glmnet_C_X4 <- unlist(lapply(glmnet_preds_X4, FUN=\(x){ coxph(Y4 ~ x)[["concordance"]]["concordance"]  }))

## SCAD performance

# C-indices
scad_C_X1 <- unlist(lapply(scad_preds_X1, FUN=\(x){ coxph(Y1 ~ x)[["concordance"]]["concordance"]  }))

# C-indices
scad_C_X2 <- unlist(lapply(scad_preds_X2, FUN=\(x){ coxph(Y2 ~ x)[["concordance"]]["concordance"]  }))

# ROC-AUCs
scad_ROCAUC_X3 <- unlist(lapply(scad_preds_X3, FUN=\(x){ auc(roc(response=Y3, predictor=x))  }))

# C-indices
scad_C_X4 <- unlist(lapply(scad_preds_X4, FUN=\(x){ coxph(Y4 ~ x)[["concordance"]]["concordance"]  }))

## APML0 performance

# C-indices
apml0_C_X1 <- unlist(lapply(apml0_preds_X1, FUN=\(x){ coxph(Y1 ~ x)[["concordance"]]["concordance"]  }))

# C-indices
apml0_C_X2 <- unlist(lapply(apml0_preds_X2, FUN=\(x){ coxph(Y2 ~ x)[["concordance"]]["concordance"]  }))

# ROC-AUCs
apml0_ROCAUC_X3 <- unlist(lapply(apml0_preds_X3, FUN=\(x){ auc(roc(response=Y3, predictor=x))  }))

# C-indices
apml0_C_X4 <- unlist(lapply(apml0_preds_X4, FUN=\(x){ coxph(Y4 ~ x)[["concordance"]]["concordance"]  }))


pdf("Fig_BigData_validations.pdf", width=12, height=5)
par(mfrow=c(1,3), las=1)
ylims <- c(0.45, 0.75)

# TCGA held out (n=101, with 18 events)
plot.new()
plot.window(xlim=c(0,50), ylim=ylims)
box(); axis(1); axis(2)
abline(h=0.5, col="grey")
abline(v=seq(from=0, to=50, by=5), col="lightgrey", lty=2)
title(xlab="Cardinality 'k' / number of non-zero variables", ylab="C-index", main="TCGA held-out data, BCR with follow-up")
# OSCAR
points(1:50, oscar_C_X4, type="l", lwd=2, col="red")
# glmnet
points(fit_glmnet$df, glmnet_C_X4, type="l", lwd=2, col="gold2")
# APML0
points(apply(bmat_apml0, MARGIN=2, FUN=\(x) { sum(!x==0) }), apml0_C_X4, type="l", lwd=2, col="blue")
# SCAD
# Omit last observation as an outlier for the highest lambda
points(apply(fit_scad$beta[,-ncol(fit_scad$beta)], MARGIN=2, FUN=\(x) { sum(!x==0) }), scad_C_X4[-ncol(fit_scad$beta)], type="l", lwd=2, col="#CC79A7")
# Greedy
points(1:50, greedy_C_X4, type="l", lwd=2, col="#009E73")

legend("bottom", col=c("red", "gold2", "blue", "#CC79A7", "#009E73", "grey"), lwd=c(2,2,2,2,2,1), legend=c("OSCAR","LASSO", "APML0", "SCAD", "Greedy FS", "Random performance"), bg="white")

# Taylor et al. independent dataset (n=131, with 27 events)
plot.new()
plot.window(xlim=c(0,50), ylim=ylims)
box(); axis(1); axis(2)
abline(h=0.5, col="grey")
abline(v=seq(from=0, to=50, by=5), col="lightgrey", lty=2)
title(xlab="Cardinality 'k' / number of non-zero variables", ylab="C-index", main="Taylor et al., BCR with follow-up")
# OSCAR
points(1:50, oscar_C_X2, type="l", lwd=2, col="red")
# glmnet
points(fit_glmnet$df, glmnet_C_X2, type="l", lwd=2, col="gold2")
# APML0
points(apply(bmat_apml0, MARGIN=2, FUN=\(x) { sum(!x==0) }), apml0_C_X2, type="l", lwd=2, col="blue")
# SCAD
# Omit last observation as an outlier for the highest lambda
points(apply(fit_scad$beta[,-ncol(fit_scad$beta)], MARGIN=2, FUN=\(x) { sum(!x==0) }), scad_C_X2[-ncol(fit_scad$beta)], type="l", lwd=2, col="#CC79A7")
# Greedy
points(1:50, greedy_C_X2, type="l", lwd=2, col="#009E73")

# Sun et al., independent dataset (n=79, of which 39 recurrent and 40 non-recurrent)
plot.new()
plot.window(xlim=c(0,50), ylim=ylims)
title(xlab="Cardinality 'k' / number of non-zero variables", ylab="ROC-AUCs", main="Sun et al., BCR binarized")
box(); axis(1); axis(2)
abline(v=seq(from=1, to=50, by=5), col="lightgrey", lty=2)
abline(h=0.5, col="grey")
# OSCAR
points(1:50, oscar_ROCAUC_X3, type="l", lwd=2, col="red")
# glmnet
points(fit_glmnet$df, glmnet_ROCAUC_X3, type="l", lwd=2, col="gold2")
# APML0
points(apply(bmat_apml0, MARGIN=2, FUN=\(x) { sum(!x==0) }), apml0_ROCAUC_X3, type="l", lwd=2, col="blue")
# SCAD
# Omit last observation as an outlier for the highest lambda
points(apply(fit_scad$beta[,-ncol(fit_scad$beta)], MARGIN=2, FUN=\(x) { sum(!x==0) }), scad_ROCAUC_X3[-ncol(fit_scad$beta)], type="l", lwd=2, col="#CC79A7")
# Greedy
points(1:50, greedy_ROCAUC_X3, type="l", lwd=2, col="#009E73")

dev.off()

save.image("validate_workspace.RData")

###
#
# Method overlap according to genes selection order
#
###

# Look at k=1:50

genes_oscar <- names(fit_oscar@kperk[[50]])

genes_greedy <- fit_greedy$selected[1:50]

genes_glmnet <- colnames(X1)[predict(fit_glmnet, s = fit_glmnet$lambda[min(which(fit_glmnet$df >= 50))], type = "nonzero")[,1]]

genes_scad <- fit_scad$beta[,min(which(apply(fit_scad$beta, MARGIN=2, FUN=\(x) sum(!x==0)) >= 50))] |>
	( \(x) { names(x)[which(!x == 0)] } )()
	
genes_apml0 <- apply(fit_apml0$Beta, MARGIN=2, FUN=\(x) { colnames(X1)[which(!x==0)] })
genes_apml0 <- genes_apml0[[min(which(unlist(lapply(genes_apml0, FUN=length)) >= 50))]]

library(ggVennDiagram)

pdf("Fig_BigData_Venn.pdf", width=8, height=8)
ggVennDiagram(
	list(
		oscar = genes_oscar,
		greedy = genes_greedy,
		glmnet = genes_glmnet,
		scad = genes_scad,
		apml0 = genes_apml0
	)
)
dev.off()



	
## Visualize some of the intersecting genes across methods
# 3 heatmaps, of each dataset
# Annotations for recurrence

g <- c("PTER", "GPRC5D", "HSPA1B", "LINC00652", "NADK", "SEC61A2", "TMC6", "LEFTY2", "SLC14A2")

X1_stat <- HeatmapAnnotation("Recurrence" = colData(mae_tcga)[c(rownames(X1), rownames(X4)),"disease_specific_recurrence_status"])
X2_stat <- HeatmapAnnotation("Recurrence" = colData(mae_taylor)[rownames(X2),"disease_specific_recurrence_status"])
X3_stat <- HeatmapAnnotation("Recurrence" = colData(mae_sun)[rownames(X3),"disease_specific_recurrence_status"])

library(ComplexHeatmap)

h1 <- ComplexHeatmap::Heatmap(t(rbind(X1[,g], X4[,g])), top_annotation = X1_stat, column_dend_reorder = order(colData(mae_tcga)[c(rownames(X1), rownames(X4)),"disease_specific_recurrence_status"]))
h2 <- ComplexHeatmap::Heatmap(t(X2[,g]), top_annotation = X2_stat, column_dend_reorder = order(colData(mae_taylor)[rownames(X2),"disease_specific_recurrence_status"]))
h3 <- ComplexHeatmap::Heatmap(t(X3[,g]), top_annotation = X3_stat, column_dend_reorder = order(colData(mae_sun)[rownames(X3),"disease_specific_recurrence_status"]))

h1 + h2 + h3

