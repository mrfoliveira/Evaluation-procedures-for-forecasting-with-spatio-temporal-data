#' Check a STARMA model for stationarity
#' 
#' @description  Checks the coefficients of a STARMA models for stationarity.
#' @details 
#' The stationarity constraints for a STAR model are the following:
#' \itemize{
#' \item STAR(2_11)
#'   \deqn{- phi_20 + abs(phi_21) <1}
#'   \deqn{abs(phi_10 + phi_11) < 1 - phi_20 - phi_21}
#'   \deqn{abs(phi_10 - phi_11) < 1 - phi_20 + phi_21}
#' \item STAR(2_10), meaning \eqn{phi_21 = 0}
#'   \deqn{phi_20 > -1}
#'   \deqn{abs(phi_10) + abs(phi_11) < 1 - phi_20}
#' \item STAR(2_00), meaning \eqn{phi_21 = 0 & phi_11 = 0}
#'   \deqn{phi_20 > -1}
#'   \deqn{abs(phi_10) < 1 - phi_20}
#' \item STAR(1_1), meaning \eqn{phi_20 = 0 & phi_21 = 0}
#'   \deqn{abs(phi_10) + abs(phi_11) < 1}
#' \item STAR(1_0), meaning \eqn{phi=20 = 0 & phi_21 = 0 & phi_11 = 0}
#'   \deqn{abs(phi_10) < 1}
#' }
#' where \code{phi_10 = model$ar[1,1]}, 
#' \code{phi_11 = model$ar[1,2]},
#' \code{phi_20 = model$ar[2,1]},
#' and \code{phi_21 = model$ar[2,2]}
#' 
#' STMA models have the same constraints for theta as the STAR contraints for phi.
#' 
#' @param model A list with components \code{ar} and \code{ma}. Each component
#' contains a matrix of the coefficients as specified in this functions' details.
#' @references \url{http://www.tandfonline.com/doi/abs/10.1080/03610918008812173}
#' @return \code{TRUE} if the model is stationary. \code{FALSE}, otherwise.
starma_stat_check <- function(model){
  
  
  # STATIONARITY CHECK
  p <- if(any(model$ar!=0)) max(which(model$ar!=0, arr.ind = T)[,1], na.rm=T) else 0
  q <- if(any(model$ma!=0)) max(which(model$ma!=0, arr.ind = T)[,1], na.rm=T) else 0
  
  assertthat::assert_that(p %in% 0:2, q %in% 0:2, p+q>0)
  order <- c(ar=p, ma=q)
  
  # (only for STAR(2_11), STAR(2_10), STAR(2_00), STAR(1_1), STAR(1_0) and same order STMAs)
  for(part in c("ar", "ma")){
    assertthat::assert_that(ncol(model[[part]])==2)
    
    phi_10 = model[[part]][1,1]
    phi_11 = model[[part]][1,2]
    phi_20 = ifelse(order[part]==2, model[[part]][2,1], 0)
    phi_21 = ifelse(order[part]==2, model[[part]][2,2], 0)
    
    # STAR(1_0)
    if(phi_20 == 0 & phi_21 == 0 & phi_11 == 0)
      if(!(abs(phi_10) < 1)) return(FALSE)
    # STAR(1_1)
    if(phi_20 == 0 & phi_21 == 0)
      if(!(abs(phi_10) + abs(phi_11) < 1)) return(FALSE)
    # STAR(2_00)
    if(phi_21 == 0 & phi_11 == 0)
      if(!(phi_20 > -1 & abs(phi_10) < 1 - phi_20)) return(FALSE)
    # STAR(2_10)
    if(phi_21 == 0)
      if(!(phi_20 > -1 & abs(phi_10) + abs(phi_11) < 1 - phi_20)) return(FALSE)
    # STAR(2_11)
    if(!(- phi_20 + abs(phi_21) <1 & abs(phi_10 + phi_11) < 1 - phi_20 - phi_21 &
         abs(phi_10 - phi_11) < 1 - phi_20 + phi_21)) return(FALSE)
  }
  return(TRUE)
}

#' Simulate spatio-temporal data using a STARMA model
#' 
#' Generate a spatio-teporal dataset according to a STARMA model.
#' @param model A list with components \code{ar} and \code{ma}. Each component
#' contains a matrix of the coefficients where the first index corresponds to the
#' row and the second index to the column.
#' @param klist A list of matrices like the ones returned by consecutive use of 
#' \code{spdep::dnearneigh} and \code{spdep::nblag} where a value higher than 0
#'  implies that the row and column locations are neighbours
#' @param n The length of the time series to be generated
#' @param rand.gen The random generator to be used. Defaults to \code{stats::rnorm}
#' @param innov A matrix of initial innovations. Defaults to 
#'  \code{matrix(rand.gen(n*ncol(klist[[1]]), ...), n, ncol(klist[[1]]))}
#' @param seed A seed to set before generating the dataset. Defaults to \code{NULL}
#' @param FUN A (possibly non-linear) function to apply to the matrix during
#'  STAR data generation.
#' @param ... Other parameters (?)
#' @return A matrix where each column contains the data for a location,
#' and each row contains the data for a time-stamp.
#' 
#' @references See example on page 3 of 
#' \url{https://cran.r-project.org/web/packages/starma/starma.pdf}
starma_sim <- function (model, klist, n, rand.gen = stats::rnorm, 
          innov = NULL, seed=NULL, FUN=function(x){x}, ...) 
{
  
  
  if(!is.null(seed)) 
    set.seed(seed)
  if(is.null(innov))
    innov <- matrix(rand.gen(n*ncol(klist[[1]]), ...), n, ncol(klist[[1]]))
  
  if (!is.list(model)) 
    stop("'model' must be list")
  if (n <= 3L) 
    stop("'n' must be strictly positive")
  assertthat::assert_that(starma_stat_check(model))
  
  p <- if(any(model$ar!=0)) max(which(model$ar!=0, arr.ind = T)[,1], na.rm=T) else 0
  q <- if(any(model$ma!=0)) max(which(model$ma!=0, arr.ind = T)[,1], na.rm=T) else 0
  
  if (length(klist)<ncol(model$ar) | length(klist)<ncol(model$ma)) 
    stop("klist length must be at length max(ncol(model$ar), ncol(model$ma))")
  
  # MATRIX CALCULATIONS
  star <- innov
  if(p>0){
    pmat <- list()
    for(p_i in 1:p){
      pmat[[p_i]] <- matrix(0, ncol(star), ncol(star))
      for(s_i in 1:ncol(model$ar)){
        if(!is.na(model$ar[p_i,s_i])) 
          pmat[[p_i]] <- pmat[[p_i]] + model$ar[p_i,s_i]*klist[[s_i]]
      }
    }
  }
  if(q>0){
    qmat <- list()
    for(q_i in 1:q){
      qmat[[q_i]] <- matrix(0, ncol(star), ncol(star))
      for(s_i in 1:ncol(model$ma)){
        if(!is.na(model$ma[q_i,s_i])) 
          qmat[[q_i]] <- qmat[[q_i]] + model$ma[q_i,s_i]*klist[[s_i]]
      }
    }
  }
  
  for (t in 3:n) {
    orig_star <- star
    # add AR matrix terms
    pstar <- matrix(0, nrow(star), ncol(star))
    if(p>0){
      for(p_i in 1:p)
        pstar[t,] <- pstar[t,] + pmat[[p_i]] %*% FUN(orig_star[t-p_i,])
    }
    # add MA matrix terms
    qstar <- matrix(0, nrow(star), ncol(star))
    if(q>0){
      for(q_i in 1:q)
        qstar[t,] <- qstar[t,] + qmat[[q_i]] %*% innov[t-q_i,] 
    }
    # add AR, MA and innovations
    star[t,] <- pstar[t,] + qstar[t,] + innov[t, ]
  }
  
  star
}


#' Exponential function
#' 
#' Calculate \code{exp(-x/C)}.A list of matrices like the ones returned by consecutive use of 
#' \code{spdep::dnearneigh} and \code{spdep::nblag} where a value higher than 0
#'  implies that the row and column locations are neighbours
#' @param x A vector
#' @param C A constant
#' @return A vector \code{exp(-x/C)}
exp_c <- function(x, C=1E4){exp(-x/C)}

#' Identity function
#' 
#' Return the identity function.
#' @param x A vector
#' @return The vector \code{x}
identity <- function(x){x}

#' Generate coefficients for a STARMA stationary process
#' 
#' @param coef_specs  A named list of coefficient specifications for the models.
#' Each slot in the list should itself contain a named list
#' with slots \code{c_10}, \code{c_11}, \code{c_20}, \code{c_21} 
#' containing either a number which will be used directly for the 
#' coefficient of that order in the AR and/or the MA components of 
#' the model, or a vector specifying an interval whithin which a 
#' coefficient will be randomly generated.
#' @param type A vector of the types of STARMA models to be used for data
#' generation. Should be \code{STAR},\code{STARMA},\code{NL_STAR} or \code{STMA}.
#' @param ndigits Number of digits to keep when rounding coefficients
#' 
#' @return A vector of stationary coefficients generated with names 
#' \code{phi_10}, \code{phi_11}, \code{phi_20}, \code{phi_21},
#' \code{theta_10}, \code{theta_11}, \code{theta_20}, \code{theta_21} and \code{FUN}.
#' 
#' @seealso \code{\link{starma_stat_check}}
generate_coef <- function(coef_specs=list(c_10=c(-2,2), c_11=c(-2,2), 
                                            c_20=c(-1,1), c_21=0),
                          type="STAR", ndigits=3){
  
  assertthat::assert_that(type %in% c("STAR","STARMA","NL_STAR","STMA"))
  assertthat::assert_that(all(names(coef_specs)==c("c_10", "c_11", "c_20", "c_21")), 
              all(unlist(lapply(coef_specs, function(x) length(x)>=1 & length(x)<3))))
  
  tofind <- 1
  while(tofind == 1){
    
    coef <- c(phi_10=0, phi_11=0, phi_20=0, phi_21=0,
              theta_10=0, theta_11=0, theta_20=0, theta_21=0,
              FUN="identity")
    
    if(type %in% c("STAR","STARMA","NL_STAR")){
      for(c in 1:length(coef_specs)){
        if(length(coef_specs[[c]])==1)
          coef[c] <- coef_specs[[c]][1]
        else
          coef[c] <- stats::runif(1, coef_specs[[c]][1], coef_specs[[c]][2])
      }
      
      if(type=="NL_STAR") 
        coef["FUN"] <- as.character(c("cos", "sin", "atan", "tanh", "exp_c")[sample(1:5, 1)])
    }
    if(type %in% c("STMA", "STARMA")){
      for(c in 1:length(coef_specs)){
        if(length(coef_specs[[c]])==1)
          coef[c+4] <- coef_specs[[c]][1]
        else
          coef[c+4] <- stats::runif(1, coef_specs[[c]][1], coef_specs[[c]][2])
      }
    }
    
    if(!is.null(round)) coef[1:8] <- round(as.numeric(coef[1:8]),digits=ndigits)
    model <- list(ar=matrix(as.numeric(coef[1:4]), 2,2, byrow=T), 
                  ma=matrix(as.numeric(coef[5:8]), 2,2, byrow=T))
    if(starma_stat_check(model)){
      tofind <- 0
    }
  }
  coef
}

#' Simulate a STARMA stationary process
#' 
#' @param Ntimes Length of time series to generate
#' @param klist A list of matrices like the ones returned by consecutive use of 
#' \code{spdep::dnearneigh} and \code{spdep::nblag} where a value higher than 0
#'  implies that the row and column locations are neighbours
#' @param coef  A named vector of stationary stationary coefficients with names 
#' \code{phi_10}, \code{phi_11}, \code{phi_20}, \code{phi_21},
#' \code{theta_10}, \code{theta_11}, \code{theta_20}, \code{theta_21} and \code{FUN}.
#' @param scale A boolean indicating whether the data should be scaled around 0.
#' @param trash Number of initial time series points to discard
#' @param seed Seed to feed to \code{starma_sim}. Default is \code{NULL}.
#' 
#' @return A list with two slots: \code{coef}, a vector of the coefficients
#'  actually used to generate the data; \code{data}, a matrix
#'  of spatio-temporal data (columns correspond to locations, rows to
#'  time-stamps).
#'  
#'  @seealso \code{\link{starma_sim}}
generate_stdata <- function(Ntimes, klist, 
                            coef, scale=FALSE, trash=100, 
                            seed=NULL){
  
  assertthat::assert_that(all(names(coef)==c("phi_10", "phi_11", "phi_20", "phi_21",
                                 "theta_10", "theta_11", "theta_20", "theta_21", "FUN")))
  
  model <- list(ar=matrix(as.numeric(coef[1:4]), 2,2, byrow=T), 
                ma=matrix(as.numeric(coef[5:8]), 2,2, byrow=T))
  assertthat::assert_that(starma_stat_check(model), 
              msg = "Coefficients should generate a stationary model.")
  
  starma <- starma_sim(model, klist, Ntimes+trash, seed=seed,
                       FUN=get(as.character(coef["FUN"])))
  starma <- starma[(trash+1):(Ntimes+trash),]

  if(scale){
    scaled_starma <- starma::stcenter(starma)
    assertthat::assert_that( (abs(max(scaled_starma)) < 1E-6 & abs(min(scaled_starma)) < 1E-6), 
                 msg = "Scaling turns data to almost all zeroes!")
    starma <- scaled_starma
  } 
  
  starma
}

