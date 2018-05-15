cat("\n\n************************************************\nEXPERIMENTS ON ARTIFICIAL DATA SETS\n************************************************\n\n")

requireNamespace("ranger", quietly=TRUE)
requireNamespace("doParallel", quietly=TRUE)

# CHANGE NUMBER OF CORES
NCORES <- 16
NUM_THREADS <- 2
# note that NCORES x NUM_THREADS will be actually used when running ranger
cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
doParallel::registerDoParallel(cores=NCORES)


DATA_PATH <-  "data/"
UTILS_PATH <- "R/"
RES_PATH <- "inst/results/artificial/"

#to_source <- list.files(UTILS_PATH, full.names=TRUE)
#for(f in to_source) source(f)

# CROSS-VALIDATION SETTINGS
XVAL_ALLOCS <- c("Trand_SPrand","Tall_SPcontig",
                 "Tall_SPrand", "Tall_SPchecker", 
                 "Tblock_SPall", "Trand_SPall", 
                 "Tblock_SPchecker", "Tblock_SPcontig",
                 "Tblock_SPrand")
XVAL_variants <- lapply(XVAL_ALLOCS, 
                        function(x) list(nfolds=16, fold.alloc.proc=x))
names(XVAL_variants) <- paste0("x", XVAL_ALLOCS)

# BUFFERED CV SETTINGS
NDXVAL_ALLOCS <- c(rep("Trand_SPrand", 3), 
                   "Tall_SPrand",
                   "Tall_SPcontig",
                   "Tblock_SPall")
NDXVAL_variants <- lapply(NDXVAL_ALLOCS, 
                          function(x) list(nfolds=16, fold.alloc.proc=x,
                                           s.buffer=1, t.buffer=3))
names(NDXVAL_variants) <- paste0("x", NDXVAL_ALLOCS)
names(NDXVAL_variants) <- paste0(names(NDXVAL_variants), c("_ST", "_T", "_S",
                                                           "_S", "_S", "_T"))
NDXVAL_variants[["xTrand_SPrand_T"]]$s.buffer <- NULL
NDXVAL_variants[["xTrand_SPrand_S"]]$t.buffer <- NULL

NDXVAL_variants <- c(NDXVAL_variants, 
                     list(xTbSr_STM=list(nfolds=16, 
                                         fold.alloc.proc="Tblock_SPrand", 
                                         s.buffer=1, t.buffer=1)))

# PREQUENTIAL SETTINGS
PRE_ALLOCS <- c("Tblock_SPall", "Tblock_SPchecker",
                "Tblock_SPcontig","Tblock_SPrand")
PRE_variants <- list()
for(win in c("growing", "sliding")){
  for(rmSP in c(TRUE, FALSE)){
    y <- lapply(PRE_ALLOCS,
                function(x) list(nfolds=16, 
                                 fold.alloc.proc=x, 
                                 window=win, 
                                 removeSP=rmSP))
    names(y) <- paste0("p", PRE_ALLOCS,"_", win)
    if(rmSP) names(y) <- paste0(names(y), "_rmSP")
    PRE_variants <- c(PRE_variants, y)
  }
}

# OTHER OOS SETTINGS
MC_variants <- list("mc47.3"=list(tr.perc=0.47,
                                         ts.perc=0.03, 
                                         nreps=16), 
                    "mc56.4"=list(tr.perc=0.56, 
                                         ts.perc=0.04, 
                                         nreps=16))

HO_variants <- list("h80.20" = list(tr.perc=0.8),
                    "h94.6" = list(tr.perc=0.94))

# EXPERIMENTAL SETTINGS
MIN_TRAIN <- 2
IN_SET_PERC <- 0.8
EVAL_FUN <- regMetrics
MODELS <- c("lm", "ranger")
EVAL_PARS <- list(eval.function=EVAL_FUN,
                  .keptTrain = TRUE, norm = TRUE)
WF_PARS <- list(lm=list(min_train = MIN_TRAIN), 
                ranger=list(min_train = MIN_TRAIN, num.threads=NUM_THREADS, verbose=FALSE))
COMPRESS <- TRUE

cat(paste("\nLoading data...\n"))
load(paste0(DATA_PATH, "lagged_data.Rdata"))

g_sizes <- names(lagged_data)
t_sizes <- names(lagged_data[[1]])
for(model in MODELS){
    for(g in 1:length(g_sizes)){
      g_size <- g_sizes[g]
      
      # get site locations
      load(paste0(DATA_PATH, "data.Rdata"))
      sites <- data[[g]][[1]][[1]]$grid$sites
      rm(data);gc()
      # get spatial distance matrix
      colnames(sites) <- c("site", "x", "y")
      sites_sf <- df2site_sf(sites, "site", "x", "y")
      s.dist <- get_spatial_dist_mat(sites_sf, site_id = "site")
      for(v in 1:length(NDXVAL_variants)){
        if(!is.null(NDXVAL_variants[[v]]$s.buffer)) 
          NDXVAL_variants[[v]][["s.dists"]] <- s.dist
      }
        
      for(t in 1:length(t_sizes)){
        t_size <- t_sizes[t]
        
        cat(paste("\nLoading data ", g_size, ", ", t_size, " ...\n"))
        load(paste0(DATA_PATH, "lagged_data.Rdata"))
        lagged_data <- list(list(lagged_data[[g_size]][[t_size]]))
        names(lagged_data) <- g_size
        names(lagged_data[[1]]) <- t_size
        gc()

        # get time distance matrix
        t.ids <- sort(unique(as.numeric(as.character(lagged_data[[g_size]][[t_size]][[1]][[1]]$time))))
        t.dist <- as.matrix(dist(t.ids))
        colnames(t.dist) <- paste0("TIME_", t.ids)
        rownames(t.dist) <- paste0("TIME_", t.ids)
        for(v in 1:length(NDXVAL_variants)){
          if(!is.null(NDXVAL_variants[[v]]$t.buffer)) 
            NDXVAL_variants[[v]][["t.dists"]] <- t.dist
        }

        ALL_variants <- c(HO_variants, MC_variants, 
                  PRE_variants, XVAL_variants, NDXVAL_variants)
        IN_ESTIMATORS <- c(rep("t_oos", length(HO_variants)),
                           rep("t_oos_mc", length(MC_variants)),
                           rep("prequential_eval", length(PRE_variants)),
                           rep("kf_xval", length(XVAL_variants)),
                           rep("nd_kf_xval", length(NDXVAL_variants)))

        if(!dir.exists(paste0(RES_PATH, model, "/"))) 
          dir.create(paste0(RES_PATH, model ,"/"), recursive=TRUE)

        cat(paste("Experiments with", model,"\n"))
        res <- run_all_experiments(models=model, 
                                      lagged_data, IN_SET_PERC, tgt~., 
                                      in_estimators=IN_ESTIMATORS, 
                                      in_est.pars=ALL_variants,
                                      out_estimator = "t_oos",
                                      out_est.pars=list(tr.perc = IN_SET_PERC),
                                      workflow = "simple_workflow", 
                                      wf.pars=WF_PARS[[model]],
                                      evaluator = "evaluate", 
                                      eval.pars=EVAL_PARS,
                                      seed=1234, site_id="site", time="time",
                                      .compress=COMPRESS, .progress=paste0(RES_PATH, model, "/res"))
        rm(res);gc();
      }
    }
}


cat("\nDone!\n\n******************\n\n")