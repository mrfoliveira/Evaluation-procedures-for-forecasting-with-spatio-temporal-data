#' Classic k-fold CV
#'
#' Fold allocation of classic k-fold CV:
#' \itemize{
#' \item shuffled time
#' \item shuffled locations
#' }
#' @param data full dataset
#' @param nfolds number of folds
#' @param time column name of time-stamp in \code{data}. 
#' Default is "time"
#' @param site_id column name of location identifier in \code{data}. 
#' Default is "site"
#' @return A list with slots:
#' \itemize{
#' \item \code{data}, possibly re-ordered
#' \item \code{f}, a vector with the fold numbers 
#' (from \code{1} to \code{nfolds}) of each row in \code{data}
#' }
#' 
#' @export
Trand_SPrand <- function(data, nfolds, time="time", site_id="site") {
  assertthat::assert_that(exists("shuffle"), exists("cv_folds"))
  
  data <- shuffle(data)
  f <- cv_folds(data, nfolds)
  
  list(data=data, f=f)
}


#' Spatially blocked CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item all time
#' \item contiguously blocked locations
#' }
#' @inheritParams Trand_SPrand
#' @inherit Trand_SPrand return
#' 
#' @export
Tall_SPcontig <- function(data, nfolds, time="time", site_id="site"){
  
  assertthat::assert_that(exists("shuffle"), exists("sp_contig"))
  
  # shuffle time
  data <- shuffle(data)
  
  # block space contiguously
  site.ids <- sort(unique(data[[site_id]]))
  sp.folds <- sp_contig(nfolds, length(site.ids))
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to fold
  f <- sp.folds[paste0("site_", data[[site_id]])]
  
  list(data=data, f=f)
}

#' Spatial CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item all time
#' \item shuffled individual locations
#' }
#' @inheritParams Trand_SPrand
#' @inherit Trand_SPrand return
#' 
#' @export
Tall_SPrand <- function(data, nfolds, time="time", site_id="site"){
  
  assertthat::assert_that(exists("shuffle"), exists("cv_folds"))
  
  # order time
  data <- data[ order(data[[time]], data[[site_id]]), ]
  
  # shuffle space
  site.ids <- shuffle(unique(data[[site_id]]))
  sp.folds <- cv_folds(site.ids, nfolds)
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to fold
  f <- sp.folds[paste0("site_", data[[site_id]])]
  
  list(data=data, f=f)
}

#' Systematic spatial CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item all time
#' \item systematically assigned (checkered) individual locations
#' }
#' @inheritParams Trand_SPrand
#' @inherit Trand_SPrand return
#' 
#' @export
Tall_SPchecker <- function(data, nfolds, time="time", site_id="site"){
  
  assertthat::assert_that(exists("sp_checker"))
  
  # shuffle time
  data <- shuffle(data)
  
  # checker space
  site.ids <- sort(unique(data[[site_id]]))
  sp.folds <- sp_checker(nfolds, length(site.ids))
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to fold
  f <- sp.folds[paste0("site_", data[[site_id]])]
  
  list(data=data, f=f)
}

#' Temporally blocked CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item blocked time
#' \item all locations
#' }
#' @inheritParams Trand_SPrand
#' @inherit Trand_SPrand return
#' 
#' @export
Tblock_SPall<- function(data, nfolds, time="time", site_id="site"){
  
  assertthat::assert_that(exists("shuffle"), exists("cv_folds"))
  
  # shuffle space
  data <- shuffle(data)
  
  # block time
  time.ids <- sort(unique(data[[time]]))
  t.folds <- cv_folds(time.ids, nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # assign rows to fold
  f <- t.folds[paste0("time_", data[[time]])]
  
  list(data=data, f=f)
}

#' Temporal CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item shuffled time
#' \item all locations
#' }
#' @inheritParams Trand_SPrand
#' @inherit Trand_SPrand return
#' 
#' @export
Trand_SPall <- function(data, nfolds, time="time", site_id="site"){
  
  assertthat::assert_that(exists("shuffle"), exists("cv_folds"))
  
  # order locations
  data <- data[ order(data[[site_id]], data[[time]]), ]
  
  # shuffle time
  time.ids <- shuffle(unique(data[[time]]))
  t.folds <- cv_folds(time.ids, nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # assign rows to fold
  f <- t.folds[paste0("time_", data[[time]])]
  
  list(data=data, f=f)
}

#' Temporal blocked and systematic spatial CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item blocked time
#' \item systematically assigned (checkered) individual locations
#' }
#' @inheritParams Trand_SPrand
#' @param t.nfolds number of folds across time. Default is sqrt(\code{nfolds})
#' @param sp.nfolds number of folds across space. Default is sqrt(\code{nfolds})
#' @return A list with slots:
#' \itemize{
#' \item \code{data}, possibly re-ordered
#' \item \code{f}, a vector with the fold identifiers of each row in \code{data}.
#' The fold identifier is composed of the concatenation of time-fold number
#' (from \code{1} to \code{t.nfolds}) and space-fold number
#' (from \code{1} to \code{sp.nfolds}), separated by "_".
#' }
#' 
#' @export
Tblock_SPchecker <- function(data, nfolds,
                               t.nfolds=round(sqrt(nfolds)), 
                               sp.nfolds=round(sqrt(nfolds)),
                             time="time", site_id="site"){
  
  assertthat::assert_that(exists("cv_folds"), exists("sp_checker"))
  
  data <- data[ order(data[[time]], data[[site_id]]), ]
  
  # block time
  time.ids <- sort(unique(data[[time]]))
  t.folds <- cv_folds(time.ids, t.nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # checker space
  site.ids <- sort(unique(data[[site_id]]))
  sp.folds <- sp_checker(sp.nfolds, length(site.ids))
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to folds
  t.f <- t.folds[paste0("time_", data[[time]])]
  sp.f <- sp.folds[paste0("site_", data[[site_id]])]
  f <- paste0(t.f, "_", sp.f)
  
  list(data=data, f=f)
}

#' Temporal blocked and randomly assigned spatial CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item blocked time
#' \item randomly assigned locations
#' }
#' @inheritParams Tblock_SPchecker
#' @inherit Tblock_SPchecker return
#' 
#' @export
Tblock_SPrand <- function(data, nfolds,
                             t.nfolds=round(sqrt(nfolds)), 
                             sp.nfolds=round(sqrt(nfolds)),
                             time="time", site_id="site"){
  
  assertthat::assert_that(exists("cv_folds"))
  
  data <- data[ order(data[[time]], data[[site_id]]), ]
  
  # block time
  time.ids <- sort(unique(data[[time]]))
  t.folds <- cv_folds(time.ids, t.nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # checker space
  site.ids <- shuffle(unique(data[[site_id]]))
  sp.folds <- cv_folds(site.ids, sp.nfolds)
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to folds
  t.f <- t.folds[paste0("time_", data[[time]])]
  sp.f <- sp.folds[paste0("site_", data[[site_id]])]
  f <- paste0(t.f, "_", sp.f)
  
  list(data=data, f=f)
}


#' Temporal blocked and contiguously-blocked spatial CV
#'
#' Fold allocation of k-fold CV using:
#' \itemize{
#' \item blocked time
#' \item contiguously blocked locations
#' }
#' @inheritParams Tblock_SPchecker
#' @inherit Tblock_SPchecker return
#' 
#' @export
Tblock_SPcontig <- function(data, nfolds,
                              t.nfolds=round(sqrt(nfolds)), 
                              sp.nfolds=round(sqrt(nfolds)),
                              time="time", site_id="site"){
  
  assertthat::assert_that(exists("cv_folds"), exists("sp_contig"))
  
  data <- data[ order(data[[time]], data[[site_id]]), ]
  
  # block time
  time.ids <- sort(unique(data[[time]]))
  t.folds <- cv_folds(time.ids, t.nfolds)
  names(t.folds) <- paste0("time_", time.ids)
  
  # checker space
  site.ids <- sort(unique(data[[site_id]]))
  sp.folds <- sp_contig(sp.nfolds, length(site.ids))
  names(sp.folds) <- paste0("site_", site.ids)
  
  # assign rows to folds
  t.f <- t.folds[paste0("time_", data[[time]])]
  sp.f <- sp.folds[paste0("site_", data[[site_id]])]
  f <- paste0(t.f, "_", sp.f)
  
  list(data=data, f=f)
}

