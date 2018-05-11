#' Calculate vector of error values
#' @describeIn calculate_vector_errors absolute error
#' 
#' @param y a vector of true values
#' @param y_hat a vector of predicted values
#' 
#' @return a vector of error values
ae <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))
  
  abs(y - y_hat)
}

#' @describeIn calculate_vector_errors squared error
se <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))
  
  (y - y_hat) ^ 2
}

#' Calculate error metrics
#' @describeIn calculate_errors mean squared error
#' @inheritParams ae
#' 
#' @return one error value
mse <- function(y, y_hat, na.rm=TRUE) mean(se(y, y_hat), na.rm = na.rm)

#' @describeIn calculate_errors mean squared error
rmse <- function(y, y_hat, na.rm=TRUE) sqrt(mse(y, y_hat, na.rm=na.rm))

#' @describeIn calculate_errors mean absolute error
mae <- function(y, y_hat, na.rm=TRUE) mean(ae(y, y_hat), na.rm = na.rm)

#' @describeIn calculate_errors normalized mean absolute error
#' 
#' @param y_train a vector of training values
nmae <- function(y, y_hat, y_train=NULL, statFUN=median, na.rm=TRUE){
  sae <- sum(ae(y, y_hat), na.rm=na.rm)
  if(!is.null(y_train)) denom <- sum(abs(y - statFUN(y_train, na.rm=na.rm)), na.rm=na.rm)
  else denom <- sum(abs(y - statFUN(y, na.rm=na.rm)), na.rm=na.rm)
  sae/denom
} 

#' @describeIn calculate_errors normalized mean squared error
nmse <- function(y, y_hat, y_train=NULL, statFUN=mean, na.rm=TRUE){
  sse <- sum(se(y, y_hat), na.rm=na.rm)
  if(!is.null(y_train)) denom <- sum((y - statFUN(y_train, na.rm=na.rm))^2, na.rm=na.rm)
  else denom <- sum((y - statFUN(y, na.rm=na.rm))^2, na.rm=na.rm)
  sse/denom
} 

#' @describeIn calculate_errors normalized root mean squared error
nrmse <- function(y, y_hat, y_train=NULL, statFUN=mean, na.rm=TRUE){
  sqrt(nmse(y, y_hat, y_train=y_train, statFUN=statFUN, na.rm=na.rm))
}

#' Calculate regression metrics
#' 
#' Calculate MAE, RMSE and utility-based regression evaluation metrics
#' @param trues a vector of true values
#' @param preds a vector of predicted values
#' @param y_train a vector of training values
#' @param norm a Boolean indicating whether to calculate normalized
#' regression metrics
#' @param aeStatFUN a function to calculate a summary of y_train for 
#' absolute error normalization. Default is median
#' @param seStatFUN a function to calculate a summary of y_train for 
#' squared error normalization. Default is mean
#' @param util a Boolean indicating whether to calculate utility-based
#' regression metrics
#' @param util.parms a named list of parameters to use for calculating utility-based
#' regression metrics. Should contain slots
#' \itemize{
#' \item \code{phi.parms} - the result of function \code{phi.control}
#' \item \code{phi.control} - if \code{phi.parms} is undefined, then \code{phi.control}
#' can be provided with a list of named arguments to feed function \code{phi.control} using
#' \code{y_train}. Default is \code{list(method = "extremes", extr.type="high")}
#' \item \code{loss.parms} - the results of function \code{loss.control}
#' \item \code{p} - Default is 0.5
#' \item \code{thr} - Relevance threshold. Default is 1
#' \item \code{beta} - Beta for F-measure. Default is 1
#' }
#' 
#' @return a named vector of calculated metrics
#' 
#' @export
regMetrics <- function(trues, preds, y_train=NULL, 
                       norm=FALSE, aeStatFUN = median, seStatFUN = mean,
                       util=FALSE, util.parms=NULL){
  
  if(length(trues)==0){
    metrics <- c(mae=NA, rmse=NA)
    if(norm) metrics <- c(metrics, nmae=NA, nmse=NA)
    if(util) metrics <- c(metrics, mu=NA, fm=NA, aucpr=NA, aucroc=NA)
  }else{
    metrics <- c(mae = mae(trues, preds),
                 rmse = rmse(trues, preds))
    
    if(norm)
      metrics <- c(metrics, nmae_tr = nmae(trues, preds, y_train, aeStatFUN),
                   nmse_tr = nmse(trues, preds, y_train, seStatFUN),
                   nrmse_tr = nrmse(trues, preds, y_train, seStatFUN),
                   nmae = nmae(trues, preds, NULL, aeStatFUN),
                   nmse = nmse(trues, preds, NULL, seStatFUN),
                   nrmse = nrmse(trues, preds, NULL, seStatFUN))
    
    if(util){
      require(uba)
      
      pP <- util.parms$phi.parms
      if(is.null(pP)){
        pPparms <- util.parms$phi.control
        if(is.null(pPparms)) pPparms <- list(method="extremes", extr.type="high")
        pP <- do.call("phi.control", c(list(y_train), pPparms))
      }
      
      lP <- util.parms$loss.parms
      if(is.null(lP)) lP <- loss.control(y_train) 
      
      p <- util.parms$p
      if(is.null(p)) p <- 0.5
      
      thr <- util.parms$thr
      if(is.null(thr)) thr <- 1
      
      beta <- util.parms$beta
      if(is.null(beta)) beta <- 1
      
      uP <- util.control(umetric="MU", p = p) ## default
      mu <- util(preds, trues, phi.parms=pP, loss.parms=lP, util.parms = uP)
      fm <- util(preds, trues, phi.parms=pP, loss.parms=lP, 
                 util.parms = list(umetric="Fm", event.thr=thr, beta=beta))
      aucpr <- util(preds, trues, phi.parms=pP, loss.parms=lP, 
                    util.parms =  util.control(umetric="AUCPR", event.thr=thr) )
      aucroc <- util(preds, trues, phi.parms=pP, loss.parms=lP, 
                     util.parms = util.control(umetric="AUCROC", event.thr=thr) ) 
      
      metrics <- c(metrics, mu=mu, fm=fm, aucpr=aucpr, aucroc=aucroc)
    }
  }

  metrics
}
