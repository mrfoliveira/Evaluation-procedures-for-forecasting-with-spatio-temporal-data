#' Time-wise holdout
#'
#' Performs one holdout experiment.
#' @param data full dataset
#' @param tr.perc percentage of data used for training.
#' Remaining will be used for testing
#' @param FUN function with arguments
#' \itemize{
#' \item \code{train} training set
#' \item \code{test} testing set
#' \item \code{time} column name of time-stamps
#' \item \code{site_id} column name of location identifiers
#' \item \code{form} a formula for model learning
#' \item \code{...} other arguments
#' }
#' @param form a formula for model learning
#' @param time column name of time-stamp in \code{data}. 
#' Default is "time"
#' @param site_id column name of location identifier in \code{data}. 
#' Default is "site_id"
#' @param .keepTrain if TRUE (default), instead of the results of 
#' \code{FUN} being directly returned, a list is created with both
#' the results and a \code{data.frame} with the time and site identifiers
#' of the observations used in the training step.
#' @param ... other arguments to FUN
#' 
#' @return The results of \code{FUN}. Usually, a data.frame
#' with location identifier \code{site_id}, time-stamp \code{time},
#' true values \code{trues} and the workflow's predictions \code{preds}.
#' 
#' @export
t_oos <- function(data, tr.perc, FUN, form,
                  time="time", site_id="site", 
                  .keepTrain = TRUE, ...){

  assertthat::assert_that(is.data.frame(data),
              time %in% colnames(data),
              site_id %in% colnames(data))
  
  # holdout percentage of data
  time.ids <- sort(unique(data[[time]]))
  tr.ids <- time.ids[1:ceiling(tr.perc*NROW(time.ids))]
  tr.inds <- which(data[[time]] %in% tr.ids)

  # split data into train and test sets
  train <- data[ tr.inds, ]
  test  <- data[ -tr.inds, ]
  
  # call function returning raw results
  hold.res <- FUN(form=form, train=train, test=test,
                  time=time, site_id=site_id, ...)
  if(.keepTrain) hold.res <- list(results = hold.res,
                            train = train[ , c(time, site_id, as.character(form[[2]])) ],
                            trainCols = colnames(train))
  hold.res
}

#' Time-wise Monte Carlo
#'
#' Performs a time-wise Monte Carlo experiment where 
#' split points are randomly chosen and a window of 
#' previous observations are used for training, with
#' a window of following observations used for testing.
#' @inheritParams t_oos
#' @param ts.perc percentage of data used for testing
#' @param nreps number of repetitions/split-points in experiment
#' @return If \code{keepTrain} is \code{TRUE}, a list where each slot
#' corresponds to one repetition or fold, containing a list with
#' slots \code{results} containing the results of \code{FUN}, and 
#' \code{train} containing a data.frame with the \code{time} and 
#' \code{site_id} identifiers of the observations used in the training
#' step. Usually, the results of \code{FUN} is a data.frame 
#' with location identifier \code{site_id}, time-stamp \code{time},
#' true values \code{trues} and the workflow's predictions
#' \code{preds}.
#' 
#' @export
t_oos_mc <- function(data, tr.perc, ts.perc, nreps, FUN, form, 
                     time="time", site_id="site", 
                     .keepTrain=TRUE, ...){
  require(assertthat)
  require(foreach)
  
  assertthat::assert_that(is.data.frame(data),
              time %in% colnames(data),
              site_id %in% colnames(data))
  
  data <- data[ order(data[[time]], data[[site_id]]), ]
  
  # get sequence of time from min to max in data
  time.ids <- unique(data[[time]])
  n <- length(time.ids)
  # get number of rows in train 
  train.size <- ceiling(tr.perc*n)
  # get number of rows in test
  test.size <- ceiling(ts.perc*n)
    
  # check that it is possible to do the specified number of repetitions
  assertthat::assert_that(n - test.size - train.size >= nreps - 1)
  
  # get range of time indices that can be split points
  selection.range <- (train.size + 1):(n - test.size + 1)
  # randomly select nreps split points
  split.points <- sort(sample(selection.range, nreps))
  
  # holdout percentage of data
  #for(i in 1:nreps){
  mc.res <- foreach(i=1:nreps) %dopar% {
    
    # get time-stamps within test set
    ts.ids <- time.ids[split.points[i]:(split.points[i] + test.size - 1)]
    # get row indices in data corresponding to test set
    ts.inds <- which(data[[time]] %in% ts.ids)
    
    # get time-stamps within training set
    tr.ids <- time.ids[(split.points[i] - train.size):(split.points[i] - 1)]
    # get row indices in data corresponding to training set
    tr.inds <- which(data[[time]] %in% tr.ids)
    
    # split data into train and test sets
    train <- data[ tr.inds, ]
    test  <- data[ ts.inds, ]
    
    # call workflow returning results
    res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)
    if(.keepTrain) res <- list(results = res,
                              train = train[ , c(time, site_id, as.character(form[[2]])) ],
                            trainCols = colnames(train))
  }
  mc.res
}

#' Cross-validation
#'
#' Performs a cross-validation experiment where folds can be 
#' allocated in different ways considering time and/or space
#' @inheritParams t_oos
#' @param nfolds number of folds for the data set to be separated into. \cr
#' If you would like to set the number of time and space folds separately, 
#' \code{nfolds} should be set to \code{NULL} and \code{t.nfolds} and
#' \code{sp.nfolds} should be fed as a list to \code{alloc.pars}
#' (only available when using \code{fold.alloc.proc} set to 
#' \code{Tblock_SPchecker}, \code{Tblock_SPcontig} or \code{Tblock_SPrand}).
#' @param fold.alloc.proc name of fold allocation function. Should be one of
#' \itemize{
#' \item \code{Trand_SPrand} -- each fold contains completely random observations.
#' The default
#' \item \code{Tall_SPcontig} - each fold includes all time and a 
#' contiguous block of space
#' \item \code{Tall_SPrand} - each fold includes all time and random 
#' locations in space
#' \item \code{Tall_SPchecker} - each fold includes all time and a 
#' set of systematically assigned (checkered) part of space
#' \item \code{Tblock_SPall} - each fold includes a block of contiguous time
#' for all locations
#' \item \code{Trand_SPall} - each fold includes random time-snapshots of
#' of all locations
#' \item \code{Tblock_SPchecker} - each fold includes a block of contiguous time
#' for a systematically assigned (checkered) part of space
#' \item \code{Tblock_SPcontig} - each fold includes a block of contiguous time
#' for a block of spatially contiguous locations
#' \item \code{Tblock_SPrand} -  each fold includes a block of contiguous time
#' for a randomly assigned part of space
#' }
#' @param alloc.pars parameters to pass onto \code{fold.alloc.proc}
#' @inherit t_oos_mc return
#' 
#' @export
kf_xval <- function(data, nfolds, FUN, form,
                    fold.alloc.proc="Trand_SPrand", alloc.pars=NULL,
                    time="time", site_id="site",
                    .keepTrain=TRUE, ...){
  require(foreach)
  
  assertthat::assert_that(is.data.frame(data),
              time %in% colnames(data),
              site_id %in% colnames(data),
              fold.alloc.proc %in% c("Trand_SPrand",
                                     "Tall_SPcontig",
                                     "Tall_SPrand",
                                     "Tall_SPchecker",
                                     "Tblock_SPall",
                                     "Trand_SPall",
                                     "Tblock_SPchecker",
                                     "Tblock_SPcontig",
                                     "Tblock_SPrand"),
              (is.null(alloc.pars) | is.list(alloc.pars)),
              .keepTrain %in% c(TRUE, FALSE))
  
  # call function that automatically allocates rows to folds
  fold_alloc <- do.call(fold.alloc.proc, c(list(data=data, nfolds=nfolds,
                                       time=time, site_id=site_id), alloc.pars))
  data <- fold_alloc$data
  folds <- fold_alloc$f
  
  assertthat::assert_that(is.vector(folds),
              if(!is.null(nfolds)) length(unique(folds)) == nfolds)
  
  cv.res <- foreach(i=unique(folds)) %dopar% {
    # each fold is used as test set once
    ts.inds <- which(folds == i)
    
    # split data into train and test sets
    train <- data[-ts.inds, ]
    test  <- data[ ts.inds, ]
    
    # call workflow returning results
    res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)
    if(.keepTrain) res <- list(results = res,
                              train = train[ , c(time, site_id, as.character(form[[2]])) ],
                            trainCols = colnames(train))
  }
  cv.res
}

#' Prequential evaluation
#'
#' Performs an evaluation procedure where training and test sets can 
#' be allocated in different ways, while always respecting the ordering 
#' provided by time (models are trained in the past and tested in the
#' relative future).
#' @inheritParams kf_xval
#' @param window type of blocked-time window ordering considered. 
#' Should be one of
#' \itemize{
#' \item \code{growing} - for each time block being tested, all previous
#'  time blocks are used for training
#' \item \code{sliding} - for each time block being tested, the immediately
#'  previous time blocks are used for training
#' }
#' @param fold.alloc.proc name of fold allocation function. Should be one of
#' \itemize{
#' \item \code{Tblock_SPall} - each fold includes a block of contiguous time
#' for all locations
#' \item \code{Tblock_SPchecker} - each fold includes a block of contiguous time
#' for a systematically assigned (checkered) part of space
#' \item \code{Tblock_SPcontig} - each fold includes a block of contiguous time
#' for a block of spatially contiguous locations
#' \item \code{Tblock_SPrand} -  each fold includes a block of contiguous time
#' for a randomly assigned part of space
#' }
#' @param removeSP argument that determines whether spatio-temporal blocks
#' including the space being used for testing should be removed from the training set.
#' Default is FALSE, meaning the information is not removed
#' @inherit t_oos_mc return
#' 
#' @export
prequential_eval <- function(data, nfolds, FUN, form,
                             window = "growing", 
                             fold.alloc.proc="Tblock_SPall", alloc.pars=NULL, 
                             removeSP = FALSE, time="time", site_id="site",
                             .keepTrain = TRUE,
                             ...){
  require(foreach)
  require(stringr)
  require(assertthat)
  
  assertthat::assert_that(is.data.frame(data),
              time %in% colnames(data),
              site_id %in% colnames(data),
              window %in% c("growing", "sliding"),
              fold.alloc.proc %in% c("Tblock_SPall",
                                     "Tblock_SPchecker",
                                     "Tblock_SPcontig",
                                     "Tblock_SPrand"),
              (is.null(alloc.pars) | is.list(alloc.pars)),
              removeSP %in% c(TRUE, FALSE),
              .keepTrain %in% c(TRUE, FALSE))
  
  # call function that automatically allocates rows to folds
  fold_alloc <- do.call(fold.alloc.proc, c(list(data=data,
                                                nfolds=nfolds,
                                                time=time, 
                                                site_id=site_id), 
                                            alloc.pars))
  data <- fold_alloc$data
  folds <- fold_alloc$f
  
  assertthat::assert_that(is.vector(folds),
              if(!is.null(nfolds)) length(unique(folds)) == nfolds)
  
  if(fold.alloc.proc != "Tblock_SPall")
    sep_folds <- apply(stringr::str_split_fixed(folds,"_", 2),
                       2, as.numeric)
  else
    sep_folds <- cbind(folds, 1)
  
    test_fold_inds <- which(sep_folds[,1]>1)
  test_folds <- sort(unique(folds[test_fold_inds]))
  
  #for(f in unique(folds)){
  pre.res <- foreach(f=test_folds) %dopar% {
      if(fold.alloc.proc != "Tblock_SPall"){
        fs <- as.numeric(stringr::str_split_fixed(f,"_", 2))
        t.f <- fs[1]
        sp.f <- fs[2]
      }else{
        t.f <- as.numeric(f)
      }
      
      # each fold is used as test set once
      ts.inds <- which(folds == f)
      
      if(window=="growing"){
        if(!removeSP | fold.alloc.proc == "Tblock_SPall") 
          tr.inds <- which(sep_folds[,1] < t.f)
        else 
          tr.inds <- which(sep_folds[,1] < t.f &
                                sep_folds[,2] != sp.f)
      }else{
        if(!removeSP | fold.alloc.proc == "Tblock_SPall") 
          tr.inds <- which(sep_folds[,1] == (t.f-1))
        else 
          tr.inds <- which(sep_folds[,1] == (t.f-1) &
                                sep_folds[,2] != sp.f)
      }
      
      # split data into train and test sets
      train <- data[ tr.inds, ]
      test  <- data[ ts.inds, ]
      
      # call workflow returning results
      res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)  
      if(.keepTrain) res <- list(results = res,
                                train = train[ , c(time, site_id, as.character(form[[2]])) ],
                            trainCols = colnames(train))
  }
  pre.res
}

#' Non-dependent cross-validation
#'
#' Performs a cross-validation experiment where folds can be 
#' allocated in different ways considering time and/or space
#' and a certain buffer around the testing set time and/or 
#' space is removed from the training set.
#' 
#' @inheritParams t_oos
#' @param nfolds number of folds for the data set to be separated into. \cr
#' If you would like to set the number of time and space folds separately, 
#' \code{nfolds} should be set to \code{NULL} and \code{t.nfolds} and
#' \code{sp.nfolds} should be fed as a list to \code{alloc.pars}
#' (only available when using \code{fold.alloc.proc} set to 
#' \code{Tblock_SPchecker}, \code{Tblock_SPcontig} or \code{Tblock_SPrand}).
#' @param t.buffer numeric value with the distance of the temporal buffer between
#' training and test sets. For each instance in the test set, instances that have 
#' a temporal distance of \code{t.buffer} or less at the same point in space are removed
#' from the training set.
#' @param s.buffer numeric value with the maximum distance of the spatial buffer between
#' training and test sets. For each instance in the test set, instances that have 
#' a spatial distance of \code{s.buffer} or less at the same point in time are removed
#' from the training set.
#' @param s.dists a matrix of the distances between the spatial IDs in \code{data}.
#' The column names and row names should be of type "SITE_<site_id>"
#' @param t.dists a matrix of the distances between the time-stamps in \code{data}.
#' The column names and row names should be of type "TIME_<time>"
#' @param fold.alloc.proc name of fold allocation function. Should be one of
#' \itemize{
#' \item \code{Trand_SPrand} -- each fold contains completely random observations.
#' The default
#' \item \code{Tall_SPcontig} - each fold includes all time and a 
#' contiguous block of space
#' \item \code{Tall_SPrand} - each fold includes all time and random 
#' locations in space
#' \item \code{Tblock_SPrand} -  each fold includes a block of contiguous time
#' for a randomly assigned part of space
#' \item \code{Tblock_SPall} - each fold includes a block of contiguous time
#' for all locations
#' }
#' @param alloc.pars parameters to pass onto \code{fold.alloc.proc}
#' @inherit t_oos_mc return
#' 
#' @export
nd_kf_xval <- function(data, nfolds, FUN, form, 
                    fold.alloc.proc="Trand_SPrand", alloc.pars=NULL,
                    t.buffer=NULL, s.buffer=NULL, s.dists=NULL, t.dists=NULL,
                    time="time", site_id="site",
                    .keepTrain=TRUE, ...){
  require(foreach)
  
  assertthat::assert_that(is.data.frame(data),
              time %in% colnames(data),
              site_id %in% colnames(data),
              fold.alloc.proc %in% c("Trand_SPrand",
                                     "Tall_SPrand",
                                     "Tall_SPcontig",
                                     "Tblock_SPall",
                                     "Tblock_SPrand"),
              (is.null(alloc.pars) | is.list(alloc.pars)),
              .keepTrain %in% c(TRUE, FALSE), msg = "Bad arguments to nd_kf_xval")
  
  assertthat::assert_that(ifelse(fold.alloc.proc %in% c("Trand_SPrand", "Tblock_SPrand"), 
                        (!is.null(t.buffer) & !is.null(t.dists)) | (!is.null(s.buffer) & !is.null(s.dists)), 
                            ifelse(grepl("Tall",fold.alloc.proc), 
                               !is.null(s.buffer) & !is.null(s.dists),
                                  !is.null(t.buffer) & !is.null(t.dists))), 
              msg = "Spatial and/or temporal buffers badly specified") 
  
  # call function that automatically allocates rows to folds
  fold_alloc <- do.call(fold.alloc.proc, c(list(data=data, nfolds=nfolds,
                                                time=time, site_id=site_id), alloc.pars))
  data <- fold_alloc$data
  folds <- fold_alloc$f
  
  assertthat::assert_that(is.vector(folds),
              if(!is.null(nfolds)) length(unique(folds)) == nfolds)
  
  cv.res <- foreach(i=unique(folds)) %dopar% {
  #for(i in folds){
    # each fold is used as test set once
    ts.inds <- which(folds == i)
    
    # split data into train and test sets
    train <- data[-ts.inds, ]
    test  <- data[ ts.inds, ]
    
    # finding which instances to remove
    if(grepl("Tall",fold.alloc.proc)){
      # spatial CV
      # remove s.buffer
      ts.ids <- unique(test[,site_id])
      buffer.ids <- unique(as.vector(unlist(sapply(ts.ids, 
                                       function(id) names(which(s.dists[paste0("SITE_", id),] <= s.buffer))))))
      buffer.ids <- substr(buffer.ids, 6, 1E6)
      cut.inds <- which(train[, site_id] %in% buffer.ids)
    }else{
      if(fold.alloc.proc=="Tblock_SPall"){
        # time-block CV
        # remove t.buffer
        ts.ids <- unique(test[,time])
        buffer.ids <- unique(as.vector(unlist(sapply(ts.ids, 
                                         function(id) names(which(t.dists[paste0("TIME_", id),] <= t.buffer))))))
        buffer.ids <- substr(buffer.ids, 6, 1E6)
        cut.inds <- which(train[, time] %in% buffer.ids)
      }else{
        if(fold.alloc.proc %in% c("Trand_SPrand", "Tblock_SPrand")){
          # standard CV
          # remove t.buffer and/or s.buffer
          ts.ids <- unique(test[,c(site_id, time)])
          cut.inds <- c()
          for(i in 1:nrow(ts.ids)){
            # ADD CUTS HERE!
            if(!is.null(t.buffer)){
              buffer.times <- names(which(t.dists[paste0("TIME_", ts.ids[i, time]),] <= t.buffer))
              buffer.times <- substr(buffer.times, 6, 1E6)
              
              cut.inds <- c(cut.inds, which(train[, time] %in% buffer.times & 
                                                train[, site_id] == ts.ids[i, site_id]))
            }
            if(!is.null(s.buffer)){
              buffer.sites <- names(which(s.dists[paste0("SITE_", ts.ids[i, site_id]),] <= s.buffer))
              buffer.sites <- substr(buffer.sites, 6, 1E6)
              cut.inds <- c(cut.inds, which(train[, site_id] %in% buffer.sites &
                                              train[, time] == ts.ids[i, time]))
            }
          }
        }
      }
    }
    
    old.train <- train
    cut.inds <- unique(cut.inds)
    # removing the instances from the training set
    if(length(cut.inds)>0){
      train <- train[-cut.inds,]
    }

    if(nrow(train)>nrow(test)){
      # call workflow returning results 
      res <- FUN(form=form, train=train, test=test, time=time, site_id=site_id, ...)
      if(.keepTrain) res <- list(results = res,
                               train = old.train[ , c(time, site_id, as.character(form[[2]])) ],
                            trainCols = colnames(train),
                            cutInds = cut.inds)  
    }else{
      res_df <- data.frame(time=character(0), station=character(0),
                           trues=numeric(0), preds=numeric(0))
      colnames(res_df)[1:2] <- c(time, site_id)
      tr_df <- data.frame(time=character(0), station=character(0),
                          value=numeric(0))
      colnames(tr_df) <- c(time, site_id, as.character(form[[2]]))
      
      if(.keepTrain) res <- list(results=res_df, train=tr_df,
                            trainCols = colnames(train))
    }

    res
      
  }

  cv.res[which(unlist(lapply(cv.res, function(x) !is.null(x))))]
}



