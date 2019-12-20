##
#
# First preliminary tests for R-code by Anni
# Main - reading data and calling functions
# 
##
##

##### Read data ######
# X = predictor matrix with observations (rows) and features (columns), dim: nrecord x nft
# Y = response matrix with observations (rows) and two columns, dim: nrecord x 2
#     - the first column is time of the event
#     - the second column is the label of the event (0=censored, 1=death) ## ONKO OIKEIN PÄIN?
# C = kit cost matrix with costs of all kits as rows, dim: nkits x 1
# K = kit matrix with all kits (rows) and feature inclusion (columns), dim: nkits x nft
#######################

library(blasso)

X <- read.table("~/blasso-Anni/data/X.txt",header=FALSE)  # Reading predictor matrix
Y <- read.table("~/blasso-Anni/data/Y.txt",header=FALSE)  # Reading response matrix
C <- read.table("~/blasso-Anni/data/C.txt",header=FALSE)  # Reading kit costs matrix
K <- read.table("~/blasso-Anni/data/K.txt",header=FALSE)  # Reading kit matrix

## Functionkutsut yms määrittelyjä

