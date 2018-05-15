cat("\n\n************************************************\nEXPERIMENTS ON REAL DATA SETS\n************************************************\n\n")

requireNamespace("ranger", quietly=TRUE)
requireNamespace("doParallel", quietly=TRUE)

# CHANGE NUMBER OF CORES
NCORES <- 9
NUM_THREADS <- 5
NUM_SPLITS <- 10
# note that NCORES x NUM_THREADS will be actually used when running ranger
cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
doParallel::registerDoParallel(cores=NCORES)

DATA_PATH <-  "data/"
UTILS_PATH <- "R/"
RES_PATH <- "inst/results/real/"
NEIB_PATH <- "R/"

if(!dir.exists(RES_PATH)){
  dir.create(RES_PATH, recursive = TRUE)
  cat(paste0("Created folder ", RES_PATH, "\n"))
} 
cat(paste0("Saving results to ", RES_PATH, "\n\n"))

#to_source <- list.files(UTILS_PATH, full.names=TRUE)
#for(f in to_source) source(f)

# CROSS-VALIDATION SETTINGS
XVAL_ALLOCS <- c("Trand_SPrand",
                 "Tall_SPrand",  
                 "Tblock_SPall",
                 "Trand_SPall", 
                 "Tblock_SPrand")
XVAL_variants <- lapply(XVAL_ALLOCS, 
                        function(x) list(nfolds=9, fold.alloc.proc=x))
names(XVAL_variants) <- paste0("x", XVAL_ALLOCS)

# BUFFERED CV SETTINGS
NDXVAL_ALLOCS <- c(rep("Trand_SPrand", 3), 
                   "Tall_SPrand",
                   "Tblock_SPall")
NDXVAL_variants <- lapply(NDXVAL_ALLOCS, 
                          function(x) list(nfolds=9, fold.alloc.proc=x,
                                           s.buffer=NA, t.buffer=NA))
names(NDXVAL_variants) <- paste0("x", NDXVAL_ALLOCS)
names(NDXVAL_variants) <- paste0(names(NDXVAL_variants), c("_ST", 
                                                           "_T", "_S", "_S", "_T"))
NDXVAL_variants[["xTrand_SPrand_T"]]$s.buffer <- NULL
NDXVAL_variants[["xTrand_SPrand_S"]]$t.buffer <- NULL

NDXVAL_variants <- c(NDXVAL_variants, list(xTbSr_STM=list(nfolds=9, fold.alloc.proc="Tblock_SPrand",
                                                     s.buffer=1, t.buffer=1)))
NDXVAL_variants <- NDXVAL_variants


# PREQUENTIAL EVALUATION SETTINGS
PRE_ALLOCS <- c("Tblock_SPall", "Tblock_SPrand")
PRE_variants <- list()
for(win in c("growing", "sliding")){
  for(rmSP in c(TRUE, FALSE)){
    y <- lapply(PRE_ALLOCS,
                function(x) list(nfolds=9, 
                                 fold.alloc.proc=x, 
                                 window=win, 
                                 removeSP=rmSP))
    names(y) <- paste0("p", PRE_ALLOCS,"_", win)
    if(rmSP) names(y) <- paste0(names(y), "_rmSP")
    PRE_variants <- c(PRE_variants, y)
  }
}

# OTHER OOS SETTINGS
MC_variants <- list("mc44.6"=list(tr.perc=0.44, 
                                   ts.perc=0.06, 
                                   nreps=9), 
                    "mc53.7"=list(tr.perc=0.53, 
                                   ts.perc=0.07, 
                                   nreps=9)) 

HO_variants <- list("h80.20" = list(tr.perc=0.8),
                    "h89.11" = list(tr.perc=0.89))

# EXPERIMENTAL SETTINGS
NORP <- 0.2
MIN_TRAIN <- 100
IN_SET_PERC <- 0.8
EVAL_FUN <- regMetrics
EVAL_PARS <- list(eval.function=EVAL_FUN, .keptTrain = TRUE, norm = TRUE)
MODELS <- c("lm", "ranger")
WF_PARS <- list(lm=list(model="lm", handleNAs="centralImput", 
                        min_train = MIN_TRAIN,
                        nORp =  NORP),
                ranger=list(model="ranger", handleNAs="centralImput", 
                            min_train = MIN_TRAIN,
                            nORp = NORP, 
                            num.threads=NUM_THREADS, verbose=FALSE))

cat("Loading data sets...\n")
load(paste0(DATA_PATH, "real_dfs.Rdata"))
dfnms <- names(data_list)

cat("Running experiments...\n\n")

res <- list()
inds_df <- list()

# for partial indicators
try(load(paste0(DATA_PATH, "inds_df.Rdata")))
for(model in MODELS){
  res[[model]] <- list()
  for(i in 1:length(data_list)){
    dfnm <- dfnms[i]
    cat(paste("\n\nExperiments with", dfnm, ":\n\n"))
   
    cat("\nCalculating spatial distance matrix...\n")
    # calculate spatial distance matrix for non-dependent X-val
    s.dist <- norm_scale(get_spatial_dist_mat(data_list[[dfnm]]$stations, site_id = "station"))
    cat("Calculating temporal distance matrix...\n")
    # calculate temporal distance matrix for non-dependent X-val
    t.dist <- norm_scale(get_time_dist_mat(data_list[[dfnm]]$df$time))
  
    ALPHA <- length(unique(data_list[[dfnm]]$df$station)) / length(unique(data_list[[dfnm]]$df$time))
    ALPHA <- ifelse(ALPHA<0.5, 0.25, 0.5)
    BETAS <- c(0.0250, 0.0375, 0.0500)
  
    # fixing NDXVAL_variants
    for(v in 1:length(NDXVAL_variants)){
      if(!is.null(NDXVAL_variants[[v]]$s.buffer)){
  		NDXVAL_variants[[v]][["s.dists"]] <- s.dist
  		if(is.na(NDXVAL_variants[[v]]$s.buffer)) NDXVAL_variants[[v]]$s.buffer <- max(BETAS)/ALPHA
      } 
    }
    for(v in 1:length(NDXVAL_variants)){
      if(!is.null(NDXVAL_variants[[v]]$t.buffer)){
  		NDXVAL_variants[[v]][["t.dists"]] <- t.dist
  		if(is.na(NDXVAL_variants[[v]]$t.buffer)) NDXVAL_variants[[v]]$t.buffer <- max(BETAS)/(1-ALPHA)
      } 
    }
  
    ALL_variants <- c(HO_variants, MC_variants, 
                    PRE_variants, XVAL_variants,
                    NDXVAL_variants)
    IN_ESTIMATORS <- c(rep("t_oos", length(HO_variants)),
                     rep("t_oos_mc", length(MC_variants)),
                     rep("prequential_eval", length(PRE_variants)),
                     rep("kf_xval", length(XVAL_variants)),
                     rep("nd_kf_xval", length(NDXVAL_variants)))
  
    if(!(dfnm %in% names(inds_df))){
      cat("Get spatio-temporal indicators...\n")
      ind_df <- get_full_indicators(data_list[[dfnm]]$df, data_list[[dfnm]]$stations,
                       k=8, var="value",
                                   betas=BETAS, alpha=ALPHA,
                                   stats = c("mean", "weighted.mean", "sd"), 
                                   ratios2add = c(TRUE,TRUE,FALSE),
                                   parallel=TRUE, nsplits=NUM_SPLITS,
                                   time_id="time", site_id="station") 
      inds_df[[dfnm]] <- list(df=ind_df, alpha=ALPHA, betas=BETAS)
  
      cat("\nSaving indicator data...\n")
      save(inds_df, file=paste0(DATA_PATH, "inds_df.Rdata"))
    }
    data_list[[dfnm]] <- NULL
    gc()
  
    # discard columns with too many NAs
    discCols <- which(sapply(
      inds_df[[dfnm]]$df[1:as.integer(IN_SET_PERC*nrow(inds_df[[dfnm]]$df)),], 
      function(y) length(which(is.na(y)))/length(y)) > NORP)
    if(length(discCols)>0){
      if(length(discCols) < 20){
        cat(paste0("Discarding columns", paste0(colnames(inds_df[[dfnm]]$df)[discCols], collapse=", "), "\n"))
        df <- as.data.frame(inds_df[[dfnm]]$df)[,-discCols]
      }else{
        cat("\n\nWARNING!\n! All columns would be discarded... None will\n\n")
        df <- as.data.frame(inds_df[[dfnm]]$df )
      }
    }else{
      df <- as.data.frame(inds_df[[dfnm]]$df )
    } 
    
    # run experiments
    cat(paste("Running experiment with", model, "...\n"))
    res[[model]][[dfnm]] <- run_one_experiment(
      df, IN_SET_PERC, as.formula("value~."), 
      in_estimators=IN_ESTIMATORS, 
      in_est.pars=ALL_variants,
      out_estimator = "t_oos",
      out_est.pars=list(tr.perc = IN_SET_PERC),
      workflow = "simple_workflow", 
      wf.pars=WF_PARS[[model]],
      evaluator = "evaluate", 
      eval.pars=EVAL_PARS,
      seed=1234, site_id="station", time="time")
    cat("Compressing results...\n")
    res[[model]][[dfnm]] <- compressRes(res[[model]][[dfnm]])
    gc()

    cat("\nSaving results...\n")
    save(res, file=paste0(RES_PATH, "real_res_R",
    	paste0(BETAS, collapse="_"),".Rdata"))
  }
}
cat("\n")

cat("\nDone!\n\n******************\n\n")