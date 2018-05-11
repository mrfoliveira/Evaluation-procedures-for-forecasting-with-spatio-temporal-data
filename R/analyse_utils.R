#' Summarize the results of one in-set/out-set experiment
#'
#' @param one_exp_res a list containing two slots: \code{out_estRes} containing a 
#' named vector of metrics estimated in out-set data, and \code{in_estRes} containing a
#' list of data frames where each column corresponds to a metric and each row to a 
#' repetition/iteration of an estimator used on in-set data
#' @param statFUN a function to summarize the evaluation metrics. Default is \code{mean}
#' @param na.rm whether to remove NAs in function \code{statFUN}
#'
#' @return A data frame with a first column containing a summary (e.g., the mean) of
#' metrics measured in the out-set data and further columns containing summaries of 
#' metrics estimated in the in-set data
#' 
#' @seealso \code{\link{run_one_experiment}}
#' 
#' @export
summarize_one_exp <- function(one_exp_res, statFUN=mean,
                              na.rm = FALSE){
  resTab <- cbind(t(one_exp_res$out_estRes$evalRes), 
                  sapply(one_exp_res$in_estRes, 
                         function(y) apply(y$evalRes, 2, statFUN, na.rm=na.rm)))
  colnames(resTab)[1] <- "real"
  resTab <- cbind(data.frame(metric=rownames(resTab), resTab))
  resTab
}

#' Summarize the results of multiple in-set/out-set experiment
#'
#' @param multi_exp_res a list containing, for each experiment, a list containing two slots: 
#' \code{out_estRes} containing a named vector of metrics estimated in out-set data, 
#' and \code{in_estRes} containing a list of data frames where each column corresponds 
#' to a metric and each row to a repetition/iteration of an estimator used on in-set data
#' @inheritParams summarize_one_exp
#'
#' @return A list of data frames of the summarized results of one experiment -- each with a 
#' first column containing a summary (e.g., the mean) of metrics measured in the out-set data
#' and further columns containing summaries of metrics estimated in the in-set data
#' 
#' @seealso \code{\link{run_multiple_experiment}}, \code{\link{summarize_one_exp}}
summarize_multiple_exp <- function(multi_exp_res, statFUN=mean,
                                   na.rm = FALSE){
  require(dplyr)
  bind_rows(lapply(multi_exp_res, summarize_one_exp), .id="lag_order")
}

#' Summarize the results of all (artificial) in-set/out-set experiments
#'
#' @param all.res a multi-level list with where the first level corresponds
#' to learning model used in the experiment, the second level contains a list
#' for each grid size of artificial data set, the third level contains a list 
#' for each time series size. Inside there is a list for each generated set,
#' the list containing a list of results for each lag embed order
#' @inheritParams summarize_one_exp
#' 
#' @return A data frame containing columns identifying the learning model,
#' grid size, time series size, type of STARMA used to generate the data,
#' order of STARMA used to generate, number of iteration of the generation process
#' with those settings, lag embed order, gold standard error (that of the out-set),
#' name of error estimator and estimated error (on the in-set)
#' 
#' @seealso \code{\link{get_best_from_multi_exp}}, \code{\link{summarize_multiple_exp}}
summarize_all_art_exps <- function(all.res, statFUN, na.rm){
  
  require(dplyr)
  sumRes <- list()
  for(model in 1:length(all.res)){
    sumRes[[model]] <- list()
    for(g.size in 1:length(all.res[[model]])){
      sumRes[[model]][[g.size]] <- list()
      for(t.size in 1:length(all.res[[model]][[g.size]])){
        sumRes[[model]][[g.size]][[t.size]] <- list()
        sumRes[[model]][[g.size]][[t.size]] <- bind_rows(lapply(all.res[[model]][[g.size]][[t.size]], function(multi.res)
          summarize_multiple_exp(multi.res, statFUN=statFUN, na.rm = na.rm)), .id="gen_model")
      }
      names(sumRes[[model]][[g.size]]) <- names(all.res[[model]][[g.size]])
    }
    names(sumRes[[model]]) <- names(all.res[[model]])
  }
  names(sumRes) <- names(all.res)
  
  sumRes <- sumRes2Tab(sumRes)
  
  sumRes
}

#' Transform a multi-level list of summarized results into a table
#'
#' @param sumRes A multi-level list of summarized results where the first level 
#' corresponds to learning model used in the experiment, the second level 
#' contains a list for each grid size of artificial data set, the third level 
#' contains a list for each time series size, and the next level contains a 
#' data frame with the results obtained in the out-set (gold-standard or 
#' "real" error) as well as estimated errors for different estimators
#' (in wide format)
#'
#' @return A data frame containing columns identifying the learning model,
#' grid size, time series size, type of STARMA used to generate the data,
#' order of STARMA used to generate, number of iteration of the generation process
#' with those settings, lag embed order, gold standard error (that of the out-set),
#' name of error estimator and estimated error (on the in-set), in long format
#' 
#' @export
sumRes2Tab <- function(sumRes){
  require(dplyr)
  require(tidyr)
  
  sumResTab <- bind_rows(lapply(sumRes, function(d)
    bind_rows(lapply(d, function(x) 
      bind_rows(x,
        .id = "t_size")),
      .id = "g_size")),
    .id = "model") %>%
    separate(gen_model, c("gen_type", "gen_order"),
             sep="\\_M\\_") %>%
    separate(gen_order, c("gen_order", "gen_it"), "\\.") %>%
    mutate(lag_order = gsub("L\\_", "", lag_order)) %>% 
    mutate_at(vars(model:metric), as.factor)
  
  sumRes_real <- sumResTab %>% select(model:real)
  sumResTab <- sumResTab %>% select(-real)
  sumRes_others <- sumResTab %>%
    gather(estimator, estimated, 9:ncol(sumResTab)) %>%
    mutate_if(is.character, as.factor)
  sumResTab <- left_join(sumRes_real, sumRes_others)
  
  sumResTab
}

#' Transform a multi-level list of summarized results into a table
#'
#' @param sumRes A multi-level list of summarized results where the first level 
#' corresponds to learning model used in the experiment, the second level 
#' contains results for each data set
#'
#' @return A data frame containing columns identifying the learning model,
#' data set "gold standard"/"real" error/ (that of the out-set), 
#' name of error estimator and estimated error (on the in-set), in long format
#' 
#' @export
realSumRes2Tab <- function(sumRes, statFUN=mean,
                              na.rm = FALSE){
  require(dplyr)
  require(tidyr)
  
  sumResTab <- bind_rows(lapply(sumRes, function(x) 
    bind_rows(lapply(x, function(y) summarize_one_exp(y, statFUN=statFUN,
                              na.rm = na.rm)), 
              .id="data")), 
    .id="model") %>%
    mutate_at(vars(model:metric), as.factor)
  
  sumRes_real <- sumResTab %>% select(model:real)
  sumResTab <- sumResTab %>% select(-real)
  sumRes_others <- sumResTab %>%
    gather(estimator, estimated, 4:ncol(sumResTab)) %>%
    mutate_if(is.character, as.factor)
  sumResTab <- left_join(sumRes_real, sumRes_others)
  
  sumResTab
}

#' Compress results from all (artificial) experiments
#'
#' @param all.res a list with multiple levels (model,
#' grid size, series size and lists of full results of multiple experiments)
#' @param rmAllRaw a boolean indicating whether the whole rawRes should
#' be removed (defaults to FALSE). If TRUE, only "train" data will be
#' removed from each set of results
#'
#' @return A multi-level list containing compressed results. Either all rawRes
#' is removed, or \code{train} is substituted by a vector of the number of instances, 
#' time and location IDs in the training set, in both \code{out_estRes} and
#' \code{in_estRes}.
#' 
#' @seealso \code{\link{summarize_all_exps}}, \code{\link{run_all_experiments}}
compressAllRes <- function(all.res, rmAllRaw=F){
  require(dplyr)
  
  for(model in 1:length(all.res)){
    for(g.size in 1:length(all.res[[model]])){
      for(t.size in 1:length(all.res[[model]][[g.size]])){
        for(df in 1:length(all.res[[model]][[g.size]][[t.size]])){
          for(l in 1:length(all.res[[model]][[g.size]][[t.size]][[df]])){
            res <- all.res[[model]][[g.size]][[t.size]][[df]][[l]]
            all.res[[model]][[g.size]][[t.size]][[df]][[l]] <- compressRes(res, rmAllRaw=rmAllRaw)
          }
        }
      }
    }
  }
  all.res
}

#' Compress results from one experiment
#'
#' @param res A list containing full results of one experiment
#' (out_estRes and in_estRes)
#' @param rmAllRaw a boolean indicating whether the whole rawRes should
#' be removed (defaults to FALSE). If TRUE, only \code{train} data will be
#' substituted by a vector of the number of instances, time and location IDs 
#' in the training set
#'
#' @return A list containing compressed results of one experiment. 
#' Either all rawRes is removed, or \code{train} substituted by a vector 
#' of the number of instances, time and location IDs in the training set
#' in both \code{out_estRes} and \code{in_estRes}
#' 
#' @seealso \code{\link{summarize_one_exps}}, \code{\link{run_one_experiment}}
compressRes <- function(res, rmAllRaw=F){
  if(rmAllRaw){
    # remove rawRes from out_est
    res$out_estRes$rawRes <- NULL
  }else{
    # remove train from rawRes in out_est
    if("train" %in% names(res$out_estRes$rawRes)){
      train <- res$out_estRes$rawRes$train
      res$out_estRes$rawRes$train <- c(times=length(unique(train[,1])), 
        stations=length(unique(train[,2])), nrows=nrow(train), 
        minTgt=min(train[,3]), meanTgt=mean(train[,3]), medTgt=median(train[,3]), maxTgt=max(train[,3]))
    }else{
      if(length(res$in_estRes[[in_est]]$rawRes)>0){
        for(f in 1:length(res$out_estRes$rawRes)){
          if("train" %in% names(res$out_estRes$rawRes[[f]])){
            train <- res$out_estRes$rawRes[[f]]$train
            res$out_estRes$rawRes[[f]]$train <- c(times=length(unique(train[,1])), 
          stations=length(unique(train[,2])), nrows=nrow(train), 
        minTgt=min(train[,3]), meanTgt=mean(train[,3]), medTgt=median(train[,3]), maxTgt=max(train[,3]))
          }
        }
      }
    } 
  }
    
  for(in_est in 1:length(res$in_estRes)){
    
    # substitute distance matrix in parameters by its summary
    if("t.dists" %in% names(res$in_estRes[[in_est]]$params))
      res$in_estRes[[in_est]]$params$t.dists <- summary(as.vector(res$in_estRes[[in_est]]$params$t.dists))
    if("s.dists" %in% names(res$in_estRes[[in_est]]$params))
      res$in_estRes[[in_est]]$params$s.dists <- summary(as.vector(res$in_estRes[[in_est]]$params$s.dists))
    
    if(rmAllRaw){
      # remove rawRes from in_est
      res$in_estRes[[in_est]]$rawRes <- NULL
    }else{
      # remove train from rawRes in in_est
      if("train" %in% names(res$in_estRes[[in_est]]$rawRes)){
        train <- res$in_estRes[[in_est]]$rawRes$train
        res$in_estRes[[in_est]]$rawRes$train <- c(times=length(unique(train[,1])), 
        stations=length(unique(train[,2])), nrows=nrow(train), 
        minTgt=min(train[,3]), meanTgt=mean(train[,3]), medTgt=median(train[,3]), maxTgt=max(train[,3]))
      }else{
        if(length(res$in_estRes[[in_est]]$rawRes)>0){
          for(f in 1:length(res$in_estRes[[in_est]]$rawRes)){
            if("train" %in% names(res$in_estRes[[in_est]]$rawRes[[f]])){
              train <- res$in_estRes[[in_est]]$rawRes[[f]]$train
              res$in_estRes[[in_est]]$rawRes[[f]]$train <- c(times=length(unique(train[,1])), 
        stations=length(unique(train[,2])), nrows=nrow(train), 
        minTgt=min(train[,3]), meanTgt=mean(train[,3]), medTgt=median(train[,3]), maxTgt=max(train[,3]))
            }
          }
        }
      }  
    }
  }
  
  res
}
  