#' Spatio-temporal data sets
#'
#' Spatio-temporal data from multiple sources used for the experimental section in 
#' "Evaluation procedures for forecasting with spatio-temporal data".
#'
#' @docType data
#'
#' @usage data(data_list)
#'
#' @format An object of class \code{list}. Each slot in list contains
#' another list with objects \code{stations} of class \code{sf}, and 
#' \code{df} of class \code{data.frame} with columns time, station, and value.
#'
#' @keywords datasets
#'
#' @references  \itemize{
#'  \item Sonja Pravilovic, Annalisa Appice, and Donato Malerba. Leveraging corre-lation across space and time to interpolate geophysical data via CoKriging.Int. J. Geogr. Inf. Sci., 32(1):191–212, 2018.
#'  \item Tomislav Hengl.GSIF: Global Soil Information Facilities, 2017.
#'  \item Caley  K  Gasch,  Tomislav  Hengl,  Benedikt  Graler,  Hanna  Meyer,  Troy  SMagney, and David J Brown.  Spatio-temporal interpolation of soil water,temperature, and electrical conductivity in 3D+ T: The Cook AgronomyFarm data set.Spat. Stat., 14:70–90, 2015.
#'  \item Edzer Pebesma. spacetime: Spatio-Temporal Data in R. J. Stat. Softw., 51(7):1–30, 2012.
#'  \item Yu Zheng, Furui Liu, and Hsun-Ping Hsieh. U-Air: When Urban Air QualityInference Meets Big Data. InProc. 19th ACM SIGKDD Int. Conf. Knowl.Discov. Data Min., KDD ’13, pages 1436–1444, New York, NY, USA, 2013.ACM.
#'  }
#'
#' @source \href{http://deohs.washington.edu/mesaair/home}{MESA}
#' \href{http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/climate-normals}{NCDC}
#' \href{http://www.tceq.state.tx.us/}{TCE}
#' \href{https://data.nal.usda.gov/project/cook-agronomy-farm}{COOK}
#' \href{http://climate.geog.udel.edu/climate/}{SAC}
#' \href{http://acm.eionet.europa.eu/databases/airbase}{RURAL}
#' \href{https://www.microsoft.com/en-us/research/publication/u-air-when-urban-air-quality-inference-meets-big-data}{BEIJ}
"data_list" 
 