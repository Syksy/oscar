##
#
# First preliminary tests for R-code by Anni
# Functions
# 
##
##


## Calling basic blasso 
## IN
## X: predictor matrix
## Y: response matrix
## K: kit matrix
## C: kit cost matrix
## problem: objective function 
##           1 = Cox's proportional hazards model with L1-penalization
##           2 = Cox's proportional hazards model with L0-penalization
## lambda: the penalization parameter in L1
## user_nk: number of nonzero elements in L0

blasso_basic <- function(X,Y,K,C,problem,lambda,user_nk){
  
  if(!is.matrix(X)) stop("Predictor matrix 'X' should be of 'matrix' class.")
  if(!is.matrix(Y)) stop("Response matrix 'Y' should be of 'matrix' class.")
  if(!is.matrix(K)) stop("Kit matrix 'K' should be of 'matrix' class.")
  if(!is.matrix(C)) stop("Kit cost matrix 'C' should be of 'matrix' class.")
  
  if(!(problem %in% c(1,2))){stop("Objective function 'problem' should be an integer 1 (for L1) or 2 (for L0).")  }
  if(user_nk <0){stop("Number of nonzero elements in L0 'user_nk' should be non-negative")}
  
  ## HUOM! Miten user_nk yläraja? Tarvitaanko? Miten se perusversio valikoi, piirteittäin vai kiteittäin.
  
  #
  nft <- ncol(X)      # Number of features
  nrecord <- nrow(X)  # Number of observations
  nkits <- nrow(K)    # Number of kits
  
  # Check dimensions
  if(nrow(Y)!=nrecord){
    stop(paste("Number of observations in the response matrix Y (",nrow(Y),") is not equal to number of observations in the predictor matrix X (",nrow(X),"). Please check both."))
  }
  if(ncol(Y)!=2){
    stop(paste("Incorrect number of columns in the response matrix Y (",ncol(Y),"). The number of columns should be 2."))
  }
  if(ncol(K)!=nft){
    stop(paste("Number of columns in kit matrix K (",ncol(K),") should be equal to the amount of features (",nft,"). Check that correct features are included."))
  }
  if(nrow(C)!=nkits){
    stop(paste("Number of kit costs (",nrow(C),") is not equal to number of kits (",nkits,"). Check that correct kits are included."))
  }
  
## Check that event is 0/1.
  if(!all(Y[,2] %in% c(0,1))){
    stop(paste("It seems that in the second column of the response matrix (Y) there is an event that is not 0 or 1. Check that the events are either 0 or 1.\n
               Additionally, check that the event is in the seconf column of response matrix."))
  }

  if(!is.double(X)) { storage.mode(X) <- 'double' }	
  if(!is.double(Y)) { storage.mode(X) <- 'double' }	
  if(!is.double(K)) { storage.mode(X) <- 'double' }	
  if(!is.double(C)) { storage.mode(X) <- 'double' }
  
  else{
    if(!is.integer(nft)) { storage.mode(nft) <- 'integer' }
    if(!is.integer(nrecord)) { storage.mode(nrecord) <- 'integer' }
    if(!is.integer(nkits)) { storage.mode(nkits) <- 'integer' }
   
    print(paste("X: nrow:", nrecord, "& ncol:", nft))
    
    ## HUOM! Matriisit/ncol/nrow varmistettava oikein päin!
    result <-.Call(c_blassocoxfit_f, as.double(Y[,2]),as.double(Y[,1]),as.double(X), as.integer(nrecord), as.integer(nft),
          as.double(K), as.double(nkits),as.double(nft),as.double(C))
  }
  if(exists(result)){
    return(result)
  }else{
    stop(paste("Result was not obtained."))
  }
}
