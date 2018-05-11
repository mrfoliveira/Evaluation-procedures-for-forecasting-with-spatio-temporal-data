#' Run one in-set/out-set error estimation experiment
#'
#' @param data a data frame
#' @param in_set_perc a fraction of the data to be used as in-set
#' @param form a learning formula
#' @param in_estimators a vector of names of estimator functions to use on the in-set data
#' @param in_est.pars a named list of the same length as \code{in_estimators} containing lists of
#' arguments to feed to each estimator apllied to the in-set data
#' @param out_estimator the name of the estimator function to use on the out-set data
#' @param out_est.pars a list containing arguments to feed to the estimator apllied to the out-set data
#' @param workflow the name of the function implementing a workflow
#' @param wf.pars a list of arguments to feed to \code{workflow}
#' @param evaluator the name of the function to calculate evaluation metrics
#' @param eval.pars a list of arguments to feed to \code{evaluator}
#' @param seed a seed to set at the start of the experiment
#' @param site_id the name of the \code{data} column containing location identifiers
#' @param time the name of the \code{data} column containing time-stamps
#'
#' @return A list with two slots: \code{out_estRes} containing the results of the out-set
#' estimators, and \code{in_estRes} containing a list of results of each estimator used on
#' in-set data
#' 
#' @export
run_one_experiment <- function(data, in_set_perc, form, 
                               in_estimators, in_est.pars, 
                               out_estimator = "t_oos", out_est.pars=list(tr.perc = in_set_perc), 
                               workflow = "simple_workflow", wf.pars=NULL, 
                               evaluator = "evaluate", eval.pars=NULL,
                               seed=1234, site_id="site", time="time"){
  
  out_est.pars$.keepTrain <- TRUE
  out_est.pars$site_id <- site_id
  out_est.pars$time <- time
  in_est.pars <- lapply(in_est.pars, function(pars){
    pars$time <- time
    pars$site_id <- site_id
    pars
  })

  # Calculate the outer estimator's results
  # (considered the real results)
  real.res <- estimates(data, form=form, estimator=out_estimator,
                        est.pars = out_est.pars, workflow = workflow, wf.pars=wf.pars, 
                        evaluator = evaluator, eval.pars=eval.pars, seed=seed) 
  real.res$params <- c(list(out_estimator=out_estimator), out_est.pars)
    
  # Get in_set from real.res
  in_set_ids <- real.res$rawRes$train
  in_set_ids <- paste0(in_set_ids[,time],"_",in_set_ids[,site_id])
  tr_dat <- data[which(paste0(data[[time]], "_", data[[site_id]]) %in% in_set_ids),]
    
  # Calculate the inner estimators' results  
  
  if(is.null(names(in_est.pars)))
    in_est_nms <- make.names(in_estimators, unique=T)
  else
    in_est_nms <- names(in_est.pars)
  perf.est <- list()
  for(j in 1:length(in_estimators)){
  	res <- estimates(tr_dat, form=form, estimator=in_estimators[j],
                     est.pars = in_est.pars[[j]], workflow = workflow, wf.pars=wf.pars, 
                     evaluator = evaluator, eval.pars=eval.pars, seed=seed)
    res$params <- in_est.pars[[j]]
    
    perf.est[[ in_est_nms[j] ]] <- res
  } 
  
  list(out_estRes = real.res, in_estRes = perf.est)
}

#'  Run multiple in-set/out-set error estimation experiments
#'
#' @param data_list a list of data frames
#' @inheritParams run_one_experiment
#'
#' @return A list containing, for each data in \code{data_list}, a list with two slots: 
#' \code{out_estRes} containing the results of the out-set estimators, and 
#' \code{in_estRes} containing a list of results of each estimator used on in-set data
#' 
#' @seealso run_one_experiment
run_multiple_experiments <- function(data_list, in_set_perc, form, 
                               in_estimators, in_est.pars, 
                               out_estimator = "t_oos", out_est.pars=list(tr.perc = in_set_perc), 
                               workflow = "simple_workflow", wf.pars=NULL, 
                               evaluator = "evaluate", eval.pars=NULL,
                               seed=1234, site_id="site", time="time"){
  
  lapply(data_list, function(dat) run_one_experiment(dat, 
                                    in_set_perc=in_set_perc, form=form, 
                                    in_estimators=in_estimators, in_est.pars=in_est.pars, 
                                    out_estimator = out_estimator, out_est.pars=out_est.pars, 
                                    workflow = workflow, wf.pars=wf.pars, 
                                    evaluator = evaluator, eval.pars=eval.pars,
                                    seed=seed, site_id=site_id, time=time))
}

#' Run multiple experiments for different grid and time series sizes
#'
#' @param models vector of names of models to use for learning
#' @param nested_data_list a nested list of data sets: top level identifies grid size,
#' second level identifies time series size, and third level contains lists of multiple data sets
#' @inheritParams run_one_experiment
#' @param .compress a Boolean indicating whether to compress results
#' @param .progress a file name to save temporary results
#' @param .verbose a Boolean indicating whether progress should be reported on during experiments
#' @param .saveMem a Boolean indicating whether partial results should be discarded
#' to save memory. Only possible if .progress is TRUE. Default is FALSE.
#' 
#' @return A nested list: top top level identifies grid size, second level identifies time
#' series size, and third level contains lists with two slots: \code{out_estRes} containing
#' the results of the out-set estimators, and \code{in_estRes} containing a list of results
#' of each estimator used on in-set data
#' 
#' @seealso \code{\link{run_multiple_experiments}}
#' 
#' @export
run_all_experiments <- function(models, nested_data_list, in_set_perc, form, 
                                in_estimators, in_est.pars, 
                                out_estimator = "t_oos", out_est.pars=list(tr.perc = in_set_perc), 
                                workflow = "simple_workflow", wf.pars=NULL, 
                                evaluator = "evaluate", eval.pars=NULL,
                                seed=1234, site_id="site", time="time",
                                .compress=FALSE, .progress=NULL, .verbose=TRUE, .saveMem=FALSE){
  
  all.res <- list()
  for(model in models){
    all.res[[model]] <- list()
    if(is.null(wf.pars)) wf.pars <- list(model=model)
    else wf.pars$model <- model
    
    for(g.size in 1:length(nested_data_list)){
      gnm <- names(nested_data_list)[g.size]
      all.res[[model]][[gnm]] <- list()
      
      for(t.size in 1:length(nested_data_list[[g.size]])){
        tnm <- names(nested_data_list[[g.size]])[t.size]
        all.res[[model]][[gnm]][[tnm]] <- list()
        
        if(.verbose) cat(paste0("\n\n", model, "; ", gnm, "; ", tnm, ": Data"))
        n_dfs <- length(nested_data_list[[g.size]][[t.size]])
        all_d_res <- foreach(d=1:n_dfs, .inorder=TRUE) %do% {
          if(.verbose) cat(paste0(" ", d))
            
          dnm <- names(nested_data_list[[g.size]][[t.size]])[d]
          data_list <- nested_data_list[[g.size]][[t.size]][[d]]
          
          res <- run_multiple_experiments(data_list, in_set_perc=in_set_perc, form=form, 
                                     in_estimators=in_estimators, in_est.pars=in_est.pars, 
                                     out_estimator = out_estimator, out_est.pars=out_est.pars, 
                                     workflow = workflow, wf.pars=wf.pars, 
                                     evaluator = evaluator, eval.pars=eval.pars,
                                     seed=seed, site_id=site_id, time=time)
          if(.compress) res <- lapply(res, compressRes) 
        }
        names(all_d_res) <- names(nested_data_list[[g.size]][[t.size]])
        all.res[[model]][[gnm]][[tnm]] <- all_d_res

        if(!is.null(.progress)){
          save(all.res, file=paste0(.progress, "_", model, "_", gnm, "_", tnm, ".Rdata"))
          if(.saveMem){
            all.res[[model]][[gnm]][[tnm]] <- NULL
            gc()
          }
        }
        if(.verbose) cat("\n\n")
      }
    }
  }
  all.res
}

