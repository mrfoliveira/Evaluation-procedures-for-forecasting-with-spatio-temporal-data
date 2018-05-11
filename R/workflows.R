#' A simple learning and prediction workflow
#' 
#' @param train a data frame for training
#' @param test a data frame for testing
#' @param time the name of the column in \code{train} and
#' \code{test} containing time-stamps
#' @param site_id the name of the column in \code{train} and
#' \code{test} containing location IDs
#' @param form a formula describing the model to learn
#' @param model the name of the algorithm to use
#' @param handleNAs string indicating how to deal with NAs.
#' If "centralImput", training observations with at least 80\%
#' of non-NA columns, will have their NAs substituted by the mean
#' value and testing observatiosn will have their NAs filled in with
#' mean value regardless.
#' @param min_train a minimum number of observations that must be
#' left to train a model. If there are not enough observations, 
#' predictions will be \code{NA}. Default is 2.
#' @param nORp a maximum number or fraction of columns with missing
#' values above which a row will be removed from train before 
#' learning the model. Only works if \code{handleNAs} was
#' set to centralImputation. Default is 0.2.
#' @param ... other parameters to feed to \code{model}
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
simple_workflow <- function(train, test, form, model="lm", 
  handleNAs=NULL, min_train=2, nORp = 0.2,
  time="time", site_id="site", ...){
  dotargs <- list(...)
  
  # get true values
  trues <- responseValues(form, test)
  
  col.inds <- which(colnames(train) %in% c(time, site_id))
  # correct default mtry if model is ranger and there is no argument given
  if(model=="ranger" & !("mtry" %in% dotargs) & is.numeric(trues))
    dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
  # pre-process NAs
  if(!is.null(handleNAs)){
    if(handleNAs=="centralImput"){
      require(DMwR2)
      idxs <- manyNAs(train, nORp = nORp)
      
      if(length(idxs)) train <- train[-idxs, ]
      if(anyNA(train)) train <- centralImputation(train)
      if(anyNA(test)) test <- centralImputation(test)
    }
  }

  if(nrow(train)>=min_train){
    # train model
    m <- do.call(model, c(list(form, train[,-col.inds]), dotargs))
    # make predictions
    preds <- if(model!="ranger") predict(m, test[,-col.inds]) else predict(m, test[,-col.inds])$predictions
    # prepare result object
    res <- data.frame(time=test[[time]], site_id=test[[site_id]],
                      trues=trues, preds=preds)
  }else{
    warning("nrow(train)<min_train", call. = FALSE)
    res <- data.frame(time=test[[time]], site_id=test[[site_id]],
                      trues=trues, preds=as.numeric(NA))
  }
  colnames(res)[1:2] <- c(time, site_id)
  res
}

#' Evalute the results of a predictive workflow
#' 
#' Calculate evaluation metrics from the raw results of a workflow
#' @param wfRes a data frame (or list of data frames) containing the results of
#' a predictive workflow with columns \code{trues} and \code{preds} containing
#' the real and predicted values, respectively
#' @param eval.function the function to be used to calculate error metrics from \code{wfRes}
#' @param .keptTrain a Boolean indicating whether \code{.keepTrain} was
#' set to \code{TRUE} in calls to estimation methods. Only useful
#' if evaluation metrics need training data.
#' @param ... parameters to pass to \code{eval.function}
#'
#' @return The results (or a list of results) of \code{eval.function} applied to 
#' the data frame (or list of data frames) in \code{wfRes}
#' 
#' @export
evaluate <- function(wfRes,
                     eval.function = get("regressionMetrics", asNamespace("performanceEstimation")),
                     .keptTrain = TRUE,
                     ...){
  
  if(!.keptTrain){
    if(!("results" %in% names(wfRes))) 
      fold.res <- t(sapply(wfRes, function(x) 
        eval.function(trues=x$results$trues, 
                      preds=x$results$preds, ...)))
    else fold.res <- t(eval.function(trues=wfRes$results$trues, 
                                     preds=wfRes$results$preds, ...)) 
  }else{
    if(!("results" %in% names(wfRes))) 
      fold.res <- t(sapply(wfRes, function(x) 
        eval.function(trues=x$results$trues, 
                      preds=x$results$preds, 
                      y_train=x$train[,3], ...)))
    else fold.res <- t(eval.function(trues=wfRes$results$trues, 
                                     preds=wfRes$results$preds, 
                                     y_train=wfRes$train[,3], ...)) 
  }
  
  fold.res
}

#' Estimate error using a chosen method
#'
#' @param data a data frame
#' @param form a formula for learning
#' @param estimator the name of an error estimator function
#' @param est.pars a named list of arguments to feed to \code{estimator}
#' @param workflow the name of the workflow to use for making predictions
#' @param wf.pars a named list of arguments to feed to \code{workflow}
#' @param evaluator the name of the function to use to calculate evaluation results
#' @param eval.pars a named list of arguments to feed to \code{evaluator}
#' @param seed a seed to set before performing estimates
#'
#' @return The results of \code{evaluator} after applying \code{estimator} to the
#' learning task
#' 
#' @export
estimates <- function(data, form, estimator="kf_xval",
                      est.pars = list(nfolds=10, 
                                      fold.alloc.proc="Trand_SPrand"), 
                      workflow = "simple_workflow", wf.pars=NULL, 
                      evaluator = "evaluate", eval.pars=NULL,
                      seed=1234){
  
  if(!is.null(seed)) set.seed(1234)
  
  res <- do.call(estimator, c(list(data=data, form=form, 
                                   FUN=get(workflow, mode="function")), 
                                est.pars, wf.pars))
  est.res <- do.call(evaluator, c(list(wfRes=res), eval.pars))
  
  list(evalRes = est.res, rawRes = res)
}


