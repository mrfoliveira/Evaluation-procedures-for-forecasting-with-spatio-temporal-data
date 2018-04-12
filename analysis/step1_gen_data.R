source("R/gen_data_utils.R")

DATA_PATH <- "data/"

## STARMA SPECIFICATIONS
Nsites <- c((8+2)^2,
            (20+2)^2)
trash <- 100
Ntimes <- c(150,300) + 3

mtypes <- c("STARMA", "STAR", "STMA", "NL_STAR")
coef_specs <- list(M_2_10 = list(c_10=c(-2,2), c_11=c(-2,2), c_20=c(-1,1), c_21=0),
                   M_2_01 = list(c_10=c(-2,0), c_11=0, c_20=c(-1,1), c_21=c(-2,1)),
                   M_2_11.1 = list(c_10=c(-1.227,0.773), c_11=c(0.733,1.277), c_20=c(-0.227,1.773), c_21=-0.733),
                   M_2_11.2 = list(c_10=c(-1.755,0.245), c_11=c(-1.755,1.755), c_20=c(-0.755,0.755), c_21=0.245))
ncoefs <- c(2,2,1,1)

cat("\nGenerating simulated data...\n")
# generate the datasets
data <- generate_multiple_datasets(Nsites, Ntimes, mtypes, 
                                   coef_specs, ncoefs, trash, 
                                   init_seed=1234, sim_seed=1234)
cat(paste("Saving to", DATA_PATH,"\n"))
save(data, file=paste0(DATA_PATH, "data.Rdata"))


## LAGGING SPECIFICATIONS
LAG_use <- c(3)
SLAGS <- list(list(c(1,1,0)))
min_time <- rep(max(LAG_use) + 1, length(LAG_use))

cat("\nEmbedding data...\n")
# lag the datasets
lagged_data <- lag_multiple_datasets(data, LAG_use, SLAGS, min_time = min_time)
cat(paste("Saving lagged data to", DATA_PATH,"\n"))
save(lagged_data, file=paste0(DATA_PATH, "lagged_data.Rdata"))
cat("\nDone!\n\n******************\n\n")