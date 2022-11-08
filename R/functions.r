

#' DOY_to_Date
#' @title Day of year to date
#' @description Takes a day of year and optionally a year, and returns a date
#' @param DOY integer day of year
#' @param year year, defaults to 2020, but having the correct year will account for leap years
#' @return a date, formatted as 'ymd'
#' @export
#'
#' @examples
DOY_to_Date <- function(DOY, year = 2020){
  start <- as.Date("01-01-2020")
  lubridate::year(start) <- year
  lubridate::yday(start) <- DOY
  return(start)
}


#' @title Day of water year to date
#' @description Takes a day of water year and optionally a year, and returns a date
#' @param DOY integer day of water year
#' @param year year, defaults to 2020, but having the correct year will account for leap years
#' @return a date, formatted as 'ymd'
#' @export
#'
#' @examples

DOWY_to_Date <- function(DOWY, year = 2020){
  start       <- as.Date("2020-09-30")
  lubridate::year(start) <- year-1
  lubridate::yday(start) <- yday(start) + DOWY
  return(start)
}

#' @title Day of snow year to date
#' @description Takes a day of snow year and optionally a year, and returns a date
#' @param DOSY integer day of snow year
#' @param year year, defaults to 2020, but having the correct year will account for leap years
#' @return a date, formatted as 'ymd'
#' @export
#'
#' @examples
DOSY_to_Date <- function(DOSY, year = 2020){
  start       <- as.Date("2020-07-31")
  lubridate::year(start) <- year-1
  lubridate::yday(start) <- yday(start) + DOSY
  return(start)
}

#' Raster_Summaries
#' @title Calculate per-layer raster summaries
#' @description Generates basic summaries for each layer of a raster stack. Calculates minimum, maximum, mean, standard deviation, and histograms, by scene, not by pixel.
#' @param Raster_Stack a raster brick containing the 12 MODIS NDVI metrics
#'
#' @return a list containing table, a dataframe of min, max, mean, and sd each layer in brick, and a list of histograms
#' @export
#'
#' @examples
Raster_Summaries <- function(Raster_Stack){
  rast_summary <- terra::summary(Raster_Stack)
  sds         <- terra::global(Raster_Stack, 'sd', na.rm = TRUE)
  means       <- terra::global(Raster_Stack, 'mean', na.rm = TRUE)
  # min_date  <- DOY_to_Date(as.integer(round(mins)))
  # max_date  <- DOY_to_Date(as.integer(round(maxes)))
  # mean_date <- DOY_to_Date(as.integer(round(means)))
  #df    <- data.frame(cbind(r_summary, sds))
  #df$metric <- rownames(df)
  histogram <- hist(Raster_Stack, plot=FALSE)
  names(histogram) <- names(Raster_Stack)
  return(list(table = rast_summary, sds = sds, histograms = histogram))
}

#' Raster_quantiles
#' @title Calculate per-layer quantiles
#' @description Generates quantiles for each layer of a raster stack, per-scene, not per-pixel.
#'
#' @param Raster_Stack a raster brick containing the 12 MODIS NDVI metrics
#' @param probs a vector of quantile values to calculate
#'
#' @return a dataframe with quantiles of each metric
#' @export
#'
#' @examples
# Raster_quantiles <- function(Raster_Stack, probs = c(0,.1,.25,.5,.75,.9,1 )){
#   quantiles <- sapply(X = c(1:11), FUN = function(X){quantile(raster(Raster_Stack, layer = X),
#                                                               probs = probs)})
#   out <- data.frame(
#     quant          = rownames(quantiles),
#     onp_SOS        = DOY_to_Date(quantiles[,1]),
#     maxp_POS       = DOY_to_Date(quantiles[,6]),
#     end_EOS        = DOY_to_Date(quantiles[,3]),
#     onv_SOS_NDVI   = quantiles[,2]  %>% round(digits = 2),
#     onv_POS_NDVI   = quantiles[,7]  %>% round(digits = 2),
#     onv_EOS_NDVI   = quantiles[,4]  %>% round(digits = 2),
#     durp_Length    = quantiles[,5]  %>% round(digits = 2),
#     tindvi         = quantiles[,11] %>% round(digits = 2),
#     ranv_rangeNDVI = quantiles[,8]  %>% round(digits = 2),
#     rtup_rate_up   = quantiles[,9]  %>% round(digits = 2),
#     rtdn_rate_dn   = quantiles[,10] %>% round(digits = 2))
#   return(out)
# }


Raster_quantiles <- function(x,  probs = c(0,.1,.25,.5,.75,.9,1 ), samp_size = 10000){
  # if(ncell(x) < samp_size){
  #   q <- global(x, quantile, probs = probs)
  # }else{
    s <- spatSample(x, size = samp_size, method = 'regular')
    qt <- sapply(s, quantile, probs = probs, na.rm=TRUE)#quantile(s, probs = probs, na.rm=TRUE)
  # }
  return(t(qt))
}

#' QC_mask
#' @title Mask layers of a raster by a specified layer of same raster
#' @description Mask layers of a raster by a specified layer of same raster
#' @param x a terra spatraster a QC layer and other layers
#' @param QC_band the layer number with QC data
#' @param QC_good_value the value of QC band that indicates good data in other layers
#'
#' @return a spatraster with all but good values masked
#'
#' @examples
QC_mask <- function(x, QC_band = 12, QC_good_value=1){
  #NDVI      <- brick(rasterfile)
  Raster_QC   <- x[[QC_band]]
  Raster_good <- mask(x, mask = Raster_QC, inverse = TRUE, maskvalues = QC_good_value)
  names(Raster_good) <- band_names
  return(Raster_good)
}

# trans_met_list <- function(input_list, start_year, end_year, band_names){
#   l <- vector(mode = 'list', length = length(band_names))
#   range <- c(1:(1+(end_year - start_year)))
#   for(band_name in band_names){
#     l[[band_name]] <- lapply(X = range, FUN = function(X){input_list[[X]][[band_name]]})
#   }
#   return(l)
# }

#' @title translate a list of rasters
#' @description takes a list of multilayer spatrasters, and returns a list of multi-layer rasters. The top-level being a list of layers, and each element of that list being a list of that layer from each input file. Useful for taking a multi-year list of multi-band rasters, and returning a multi-band list of multi-year rasters.
#' @param input_list a list of multi-layer spatrasters with the same layer structure
#' @param band_names a vector of the names of the layers shared by each input file, will become the names of each element of the top-level list
#' @return a list of lists
Trans_rast_list <- function(input_list,  layer_names){
  l <- vector(mode = 'list')
  for(layer_name in layer_names){
    l[[layer_name]] <- lapply(X = input_list, FUN = function(X){X[[layer_name]]}) %>% terra::rast()
  }
  return(l)
}


