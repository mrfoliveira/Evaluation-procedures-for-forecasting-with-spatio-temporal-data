#' Feature scaling
#' 
#' Normalize values to be within the range between [0,1].
#' @param x a vector of values
#' @return a scaled vector
norm_scale <- function(x){
  if(min(x)!=max(x)) (x - min (x)) / ( max(x) - min(x) )
  else x
} 

#' Get response values of a dataset from a formula
#' 
#' @param formula learning formula
#' @param data data set to get the target values from
#' @param na what action to perform if NAs are present. Default is na.fail
#' @return A vector of the target values.
responseValues <- function (formula, data, na = NULL) 
  stats::model.response(stats::model.frame(formula, data, na.action = na))

#' Shuffle values/rows
#' 
#' Shuffle the values or rows of a vector or data frame
#' @param x a vector or data frame
#' @return a vector or data frame
shuffle <- function(x){
  if(is.null(dim(x))) x[sample(NROW(x))]
  else x[sample(NROW(x)),]
}

#' Cut into folds
#' 
#' Assigns rows of a data frame into folds for cross-validation.
#' @param x a data.frame
#' @param nfolds number of folds
#' @return a vector with the fold assignment of each row
cv_folds <- function(x, nfolds) {
  cut(seq_len(NROW(x)), breaks = nfolds, labels = FALSE)
}

#' Assign the locations of a regular grid to folds
#' following a checkered pattern
#' 
#' Systematically assigns the locations of a data frame
#' into folds which are checkered across space for 
#' cross-validation. Assumes the sites are
#' sorted (e.g., left to right, bottom to top).
#' @param nfolds number of folds to divide the space into
#' @param nsites number of locations in the regular grid 
#' @param grid.h height of the grid (in number of sites). 
#' Default is \code{sqrt(nfolds)}
#' @param grid.w width of the grid (in number of sites). 
#' Default is \code{sqrt(nfolds)}
#' @return a vector with the fold assignment of each location
sp_checker <- function(nfolds, nsites, grid.h=sqrt(nsites), grid.w=sqrt(nsites)){
  # require(wavethresh)
  
  if(nfolds<grid.h){
    nreps <- floor(grid.h/nfolds)
    remainder <- grid.h %% nfolds
    seq2rep <- c(rep(1:nfolds, nreps), 1:remainder)[1:grid.h]
  }else{
    nreps <- ceiling(nsites/nfolds)
    remainder <- 0
    seq2rep <- 1:nfolds
  }
  as.vector(sapply(0:(grid.h-1), function(i) wavethresh::guyrot(seq2rep, n=i)))[1:nsites]
}

#' Assign the locations of a regular grid to folds
#' in contiguous square blocks.
#' 
#' Assigns the locations of a data frame into contiguous 
#' blocks folds for cross-validation. Assumes the sites are
#'  sorted (e.g., left to right, bottom to top).
#' WARNING: Works well for perfect squares that can be 
#' divided into \code{nfolds} perfect squares ONLY.
#' @inheritParams sp_checker
#' @return a vector with the fold assignment of each location
sp_contig <- function(nfolds, nsites, grid.h=sqrt(nsites), grid.w=sqrt(nsites)){
  mat <- matrix(0, grid.h, grid.w)
  nhood_side <- floor(sqrt((grid.w*grid.h)/nfolds))
  folds <- paste(ceiling(col(mat)/nhood_side), ceiling(row(mat)/nhood_side), sep="-")
  folds <- factor(folds, labels=1:length(unique(folds)))
  blocks <- matrix(folds, ncol=ncol(mat)) 
  as.vector(blocks) 
}

