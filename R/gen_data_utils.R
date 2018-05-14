source("R/starma_utils.R")

#' Get directional neighbours from a regular grid
#' 
#' Given a matrix detailing the position of location in a 
#' regular grid, returns a matrix specifying the IDs of 
#' immediate neighbours to the right, left, top and bottom 
#' of each location.
#' @param sites a matrix or data.frame with site ID equal to row number and 
#' two columns, corresponding to position in x and y
#' @param klist a matrix like the ones returned by consecutive use of 
#' spdep::dnearneigh and spdep::nblag where a value higher than 0 implies 
#' that the row and column locations are neighbours
#' @return A matrix with \code{nrow(sites)} rows and four columns containing 
#' the ID of a location's neighbour in each cardinal direction 
#' (top, right, bottom and left).
grid_neibs_ord1 <- function(sites, klist){
  # get directional neighbour IDs (order1)
  neibs_order1 <- as.data.frame(matrix(NA, nrow(sites), 4))
  colnames(neibs_order1) <- c("tneib", "rneib", "bneib", "lneib")
  if(is.data.frame(sites)) sites <- as.matrix(sites[,c("x", "y")])
  
  for(s in 1:nrow(sites)){
    site <- sites[s,]
    neibs <- which(klist[s,]>0)
    neib_sites <- sites[neibs,]
    
    tneib <- neibs[which(neib_sites[,1]==site[1] & neib_sites[,2]>site[2])]
    if(length(tneib)) neibs_order1[s,"tneib"] <- tneib
    rneib <- neibs[which(neib_sites[,1]>site[1] & neib_sites[,2]==site[2])]
    if(length(rneib)) neibs_order1[s,"rneib"] <- rneib
    bneib <- neibs[which(neib_sites[,1]==site[1] & neib_sites[,2]<site[2])]
    if(length(bneib)) neibs_order1[s,"bneib"] <- bneib
    lneib <- neibs[which(neib_sites[,1]<site[1] & neib_sites[,2]==site[2])]
    if(length(lneib)) neibs_order1[s,"lneib"] <- lneib
  }
  neibs_order1
}

#' Create a spatio-temporal embed of data
#' 
#' Given spatio-temporal data, transform it into a data frame
#' where each row has its own value as target variable and,
#' as predictors, lagged values recorded in the past at 
#' its own location and/or in the past at neighbouring 
#' locations.
#' @param data a matrix of spatio-temporal data where each
#'  row contains the variable information for a time point
#'  and each column corresponds to a certain location
#' @param neibs a data frame with information about the
#' neighbours, like the one returned by \code{grid_neibs_ord1}
#' @param p a number for the order of the temporal lag for 
#' values at a specific location
#' @param slags a vector with the length of the order of
#' temporal lags for values at immediately neighbouring
#' locations. A value of 0 means no neighbours are used at 
#' that temporal lag, a value of 1 means that the values
#' of immediate neighbours are included for that temporal lag.
#' @return a data frame with columns time, site, target, 
#' as many self_lag columns as \code{p}, and as many
#' top, bottom, right and left neighbour value lags as the 
#' sum of \code{slags}.
#' @seealso \code{\link{grid_neibs_ord1}}
st_lag_neib_ord1 <- function(data, neibs, p, slags){
  ncols <- sum(slags*4) + p + 1
  colnms <- c()
  for(i in 1:p){
    if(slags[i]<=1) colnms <- c(colnms, paste0("self_lag", i))
    if(slags[i]==1) colnms <- c(colnms, paste0(c("t_lag", "r_lag", "b_lag", "l_lag"), i))
  }
  newdata <- as.data.frame(matrix(NA, ncol=ncols+2, nrow=ncol(data)*(nrow(data)-2)))
  colnames(newdata) <- c("site", "time", "tgt", colnms)
  
  for(s in 1:ncol(data)){
    srows <- ((s-1)*(nrow(data)-p)+1):(s*(nrow(data)-p))
    newdata[srows,"site"] <- s
    newdata[srows,"time"] <- (p+1):nrow(data)
    newdata[srows, "tgt"] <- data[(p+1):nrow(data),s]
    
    for(i in 1:p){
      if(slags[i] <= 1){
        lag_line <- Hmisc::Lag(data[,s], shift=i)
        newdata[srows, paste0("self_lag", i)] <- lag_line[(p+1):length(lag_line)]
      }
      if(slags[i]==1){
        for(dir in c("t", "r", "b", "l")){
          neib <- neibs[s, paste0(dir, "neib")]
          lag_line <- Hmisc::Lag(data[,neib], shift=i)
          newdata[srows, paste0(dir, "_lag", i)] <- lag_line[(p+1):length(lag_line)]
        }
      }
    }
  }
  newdata
}

#' Create a regular grid using spdep package
#' 
#' Create a regular grid with a certain number of locations,
#' returning each location's position and neighbour matrices
#' of order 0 and 1.
#' @param Nsites Number of locations in the grid
#' @param grid.h Height of the grid (in number of locations).
#' Defaults to \code{sqrt(Nsites)}
#' @param grid.w Width of the grid (in number of locations).
#' Defaults to \code{sqrt(Nsites)}
#' @return A list with slots sites -- a data frame with
#' columns id, x and y --, and klist -- a list of neighbour
#' matrices of order 0 and 1 where the values are calculated
#' using functions from \code{spdep} package. 
#' @seealso \code{\link{dnearneigh}}, \code{\link{nblag}},
#' \code{\link[spdep]{nb2mat}}
generate_grid <- function(Nsites, grid.h=ceiling(sqrt(Nsites)),
                          grid.w=ceiling(sqrt(Nsites))){
  
  # require(spdep)
  sites <- matrix(0, grid.h*grid.w, 2)
  for (i in 1:grid.h) {
    for (j in 1:grid.w)
      sites[(i-1)*grid.w + j, ] <- c(i, j) - .5
  }
  
  # Create a uniform first order neighbourhood
  knb <- spdep::dnearneigh(sites, 0, 1)
  
  # Lag the neighbourhood to create other order matrices
  knb <- spdep::nblag(knb, 2)
  
  klist <- list(order0=diag(Nsites),
                order1=spdep::nb2mat(knb[[1]]))
  
  sites <- as.data.frame(cbind(1:NROW(sites), sites))
  colnames(sites) <- c("id", "x", "y")
  
  list(sites=sites, klist=klist)
}

#' Generate an artificial dataset using STARMA on a regular grid
#' 
#' This function generates a regular grid and simulates a time series
#' for each location according to a STARMA model specification.
#' @inheritParams generate_grid
#' @param Ntimes Number of points in each time series to generate
#' @param coef_spec Coefficient specifications. Should be a named list
#' with slots \code{c_10}, \code{c_11}, \code{c_20}, \code{c_21} 
#' containing either a number which will be used directly for the 
#' coefficient of that order in the AR and/or the MA components of 
#' the model, or a vector specifying an interval whithin which a 
#' coefficient will be randomly generated.
#' @param gen_coef A Boolean indicating whether \code{coef_spec} should be used
#' to generate coefficients (TRUE: default) or whether \code{coef_spec} already
#' contains a proper list of coefficients
#' @param mtype Type of STARMA model (should only be provided if \code{gen_coef} 
#' is \code{TRUE}. One of \code{STAR},\code{STARMA},\code{NL_STAR} or \code{STMA}.
#' @param trash A number of initial time series points to be discarded
#' at each location.
#' @param ndigits Number of digits to keep when rounding coefficients.
#' @param seed Seed to feed to \code{generate_stdata}. Default is \code{NULL}.
#' 
#' @return A list with three slots: \code{data}, a matrix of spatio-temporal data 
#' (columns correspond to locations, rows to time-stamps); \code{coef}, a vector of 
#' the coefficients actually used to generate the data; and \code{grid}, a list
#' containing the positions of location and the neighbour matrices.
#'  
#'  @seealso \code{\link{generate_stdata}}, \code{\link{generate_coef}}, \code{\link{generate_grid}}
generate_one_dataset <- function(Nsites, Ntimes, coef_spec, gen_coef = TRUE, mtype=NULL, 
                                 trash=0, ndigits=3, seed=NULL,
                                 grid.h=ceiling(sqrt(Nsites)), grid.w=ceiling(sqrt(Nsites))){
  
  # Get grid sites and klist
  grid <- generate_grid(Nsites, grid.h = grid.h, grid.w = grid.w)
  
  # generate starma data, getting list of coefficients used and data
  if(gen_coef) coef_spec <- generate_coef(coef_specs=coef_spec, type=mtype, ndigits=ndigits)
  dat <- generate_stdata(Ntimes, grid$klist, coef=coef_spec, trash = trash, seed = seed)
  
  list(data = dat, coef = coef_spec, grid=grid)
}

#' Generate multiple artificial datasets using STARMA on regular grids
#' 
#' This function generates multiple regular grids and simulates time series
#' for each location according to one or more STARMA model specifications.
#' @inheritParams generate_one_dataset
#' @param Nsites A vector containing the number of locations in the regular
#' grids to be generated
#' @param Ntimes A vector containing the number of time series points to be
#' generated for each location
#' @param mtypes A vector of the types of STARMA models to be used for data
#' generation. Should be \code{STAR},\code{STARMA},\code{NL_STAR} or \code{STMA}.
#' @param coef_specs A named list of coefficient specifications for the models.
#' Each slot in the list should itself contain a named list
#' with slots \code{c_10}, \code{c_11}, \code{c_20}, \code{c_21} 
#' containing either a number which will be used directly for the 
#' coefficient of that order in the AR and/or the MA components of 
#' the model, or a vector specifying an interval whithin which a 
#' coefficient will be randomly generated. Coefficients are only randomly generated
#' for the first pair of \code{Nsites} and \code{Ntimes}, and then are re-used
#' for different grid and time series sizes.
#' @param ncoefs A vector with the number of data sets that should be generated
#'  according to each of the specifications in the \code{coef_specs} list.
#' @param trash A number of initial time series points to be discarded at 
#' each location.
#' @param grid.h Height of the grid (in number of locations).
#' Defaults to \code{sqrt(Nsites)}.
#' @param grid.w Width of the grid (in number of locations).
#' Defaults to \code{sqrt(Nsites)}.
#' @param init_seed Seed to set at the begining (before coefficient generation).
#' Default is \code{1234}.
#' @param mid_seed  Seed to set between grid/time series size change.
#' Default is \code{NULL}.
#' @param sim_seed Seed to feed to \code{generate_one_dataset}.
#' Default is \code{NULL}.
#' 
#' 
#' @return A list containing a list for each size in \code{Nsites}. Each
#' of the sublists contains a list of datasets generated by 
#' \code{generate_one_dataset} for each time series size specified in 
#' \code{Ntimes}. Each of these sublists will have 
#' \code{sum(length(mtypes) x ncoefs)} datasets of  each grid and
#'  time series size. In total, \code{length(Nsites) x length(Ntimes) 
#'  x sum(length(mtypes) x ncoefs)} spatio-temporal datasets are generated.
#' @seealso \code{\link{generate_one_dataset}}
#' 
#' @export
generate_multiple_datasets <- function(Nsites, Ntimes, mtypes, coef_specs, ncoefs, trash=0,
                                       grid.h = ceiling(sqrt(Nsites)), 
                                       grid.w = ceiling(sqrt(Nsites)),
                                       init_seed = 1234, mid_seed = NULL, sim_seed = NULL){
  
  if(!is.null(init_seed))
    set.seed(init_seed)
  
  coefs <- list()
  for(c in 1:length(coef_specs)){
    for(mtype in 1:length(mtypes)){
      for(i in 1:ncoefs[c]){
        
        nm <- paste0(mtypes[mtype], "_", names(coef_specs)[c])
        if(ncoefs[c] > 1) nm <- paste0(nm, ".", i)
        coefs[[nm]] <- generate_coef(coef_specs=coef_specs[[c]], type=mtypes[mtype])
      }
    }
  }
  
  data <- list()
  for(g.size in 1:length(Nsites)){
    data[[g.size]] <- list()
    for(t.size in 1:length(Ntimes)){
      data[[g.size]][[t.size]] <- list()
      if(!is.null(mid_seed)) set.seed(mid_seed)
      for(c in 1:length(coefs)){
        dat <- generate_one_dataset(Nsites = Nsites[g.size], Ntimes = Ntimes[t.size], 
                                    coef_spec = coefs[[c]], gen_coef = FALSE,
                                    trash=trash, seed = sim_seed,
                                    grid.h = grid.h[g.size], grid.w = grid.w[g.size])
        data[[g.size]][[t.size]][[c]] <- dat
      }
      names(data[[g.size]][[t.size]]) <- names(coefs)
    }
    names(data[[g.size]]) <- paste0("TIME_SZ_", unlist(lapply(data[[g.size]], 
                                                              function(y) 
                                                                nrow(y[[1]]$data))))
  }
  names(data) <- paste0("GRID_SZ_", unlist(lapply(data, 
                                                  function(y) 
                                                    ncol(y[[1]][[1]]$data))))
  data
}
  
#' Create a spatio-temporal embed of multiple datasets
#' 
#' Transform multiple spatio-temporal data sets in a list of lists
#' into a list of lists containing data frames that have a
#' spatio-temporal embed. Only data for "interior" points in the grid 
#' is kept. That is, points that do not have four immediate neighbours 
#' are filtered out.
#' @param data_list A multi-level list as the ones generated by
#' \code{generate_multiple_datasets}. The first level corresponds to 
#' a certain grid size; the second level to a certain time series size
#' and the third level contains multiple lists containing datasets (that might 
#' have been generated from different STARMA specifications).
#' @param LAG_use A vector of the temporal lag orders to use.
#' @param SLAGS A list containing vectors of length of the corresponing 
#' temporal lag size. Each of the vectors contains values of 0 meaning 
#' no neighbour columns are used at that temporal lag, and/or values of 1 
#' meaning that the values of immediate neighbours are included for that
#' temporal lag.
#' @param min_time A vector of the minimum time-stamps to be included in
#'  the datasets. Defaults to \code{rep(max(LAG_use), length(LAG_use)}.
#' @return A list of lists similar to data_list, but with the data sets
#' having spatio-temporal embeds.
#' @seealso \code{\link{generate_multiple_datasets}}, 
#' \code{\link{st_lag_neib_ord1}}
#' 
#' @export
lag_multiple_datasets <- function(data_list, LAG_use, SLAGS, 
                                  min_time=rep(max(LAG_use), length(LAG_use))){
  lagged_data <- list()
  assertthat::assert_that(length(LAG_use) == length(min_time))
  
  for(g.size in 1:length(data_list)){
    lagged_data[[g.size]] <- list()
    for(t.size in 1:length(data_list[[g.size]])){
      lagged_data[[g.size]][[t.size]] <- list()
      for(df in 1:length(data_list[[g.size]][[t.size]])){
        
        dat <- data_list[[g.size]][[t.size]][[df]]
        
        lag_dat_list <- list()
        for(l in 1:length(LAG_use)){
          lag <- LAG_use[l]
          slags <- SLAGS[[l]]
          for(s in 1:length(slags)){
            nme <- paste0("L_", lag, "_", paste0(slags[[s]], collapse=""))
            
            neibs_order1 <- grid_neibs_ord1(dat$grid$sites, dat$grid$klist[[2]])
            # check which sites have complete neighbourhoods
            int_points <- dat$grid$sites[stats::complete.cases(neibs_order1),]
            
            lag_dat <- st_lag_neib_ord1(dat$data, neibs_order1, lag, slags[[s]])
            lag_dat <- lag_dat[which(lag_dat$site %in% int_points$id),]
            lag_dat <- lag_dat[which(lag_dat$time>=min_time[l]),]
            lag_dat_list[[nme]] <- lag_dat
          }
        }
        df_name <- names(data_list[[g.size]][[t.size]])[df]
        lagged_data[[g.size]][[t.size]][[df_name]] <- lag_dat_list
      }
    }
    names(lagged_data[[g.size]]) <- paste0("TIME_SZ_", 
                                           unlist(lapply(lagged_data[[g.size]], 
                                                         function(y) 
                                                           length(unique(y[[1]][[1]]$time)))))
  }
  names(lagged_data) <- paste0("GRID_SZ_", 
                               unlist(lapply(lagged_data, 
                                                         function(y) 
                                                           length(unique(y[[1]][[1]][[1]]$site)))))
  lagged_data
}



