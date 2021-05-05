

#' DOY_to_Date
#' Takes a day of year and optionally a year, and returns a date
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


#' DOWY_to_Date
#' Takes a day of water year and optionally a year, and returns a date
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

#' DOSY_to_Date
#' Takes a day of snow year and optionally a year, and returns a date
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
#' Generates basic summaries for each layer of a raster stack. Calculates minimum, maximum, mean, standard deviation, and histograms, by scene, not by pixel.
#' @param Raster_Stack a raster brick containing the 12 MODIS NDVI metrics
#'
#' @return a list containing table, a dataframe of min, max, mean, and sd each layer in brick, and a list of histograms
#' @export
#'
#' @examples
Raster_Summaries <- function(Raster_Stack){
  mins  <- summary(Raster_Stack, 'min')
  maxes <- cellStats(Raster_Stack, 'max')
  means <- cellStats(Raster_Stack, 'mean')
  sds   <- cellStats(Raster_Stack, 'sd')
  # min_date  <- DOY_to_Date(as.integer(round(mins)))
  # max_date  <- DOY_to_Date(as.integer(round(maxes)))
  # mean_date <- DOY_to_Date(as.integer(round(means)))
  df    <- data.frame(mins, maxes, means, sds)
  df$metric <- rownames(df)
  histogram <- hist(Raster_Stack, plot=FALSE)
  return(list(table = df, histograms = histogram))
}

#' Raster_quantiles
#' Generates quantiles for each layer of a raster stack, per-scene, not per-pixel.
#'
#' @param Raster_Stack a raster brick containing the 12 MODIS NDVI metrics
#' @param probs a vector of quantile values to calculate
#'
#' @return a dataframe with quantiles of each metric
#' @export
#'
#' @examples
Raster_quantiles <- function(Raster_Stack, probs = c(0,.1,.25,.5,.75,.9,1 )){
  quantiles <- sapply(X = c(1:11), FUN = function(X){quantile(raster(Raster_Stack, layer = X),
                                                              probs = probs)})
  out <- data.frame(
    quant          = rownames(quantiles),
    onp_SOS        = DOY_to_Date(quantiles[,1]),
    maxp_POS       = DOY_to_Date(quantiles[,6]),
    end_EOS        = DOY_to_Date(quantiles[,3]),
    onv_SOS_NDVI   = quantiles[,2]  %>% round(digits = 2),
    onv_POS_NDVI   = quantiles[,7]  %>% round(digits = 2),
    onv_EOS_NDVI   = quantiles[,4]  %>% round(digits = 2),
    durp_Length    = quantiles[,5]  %>% round(digits = 2),
    tindvi         = quantiles[,11] %>% round(digits = 2),
    ranv_rangeNDVI = quantiles[,8]  %>% round(digits = 2),
    rtup_rate_up   = quantiles[,9]  %>% round(digits = 2),
    rtdn_rate_dn   = quantiles[,10] %>% round(digits = 2))
  return(out)
}

QC_mask <- function(rasterfile, QC_band = 12, QC_good_value=1){
  #NDVI      <- brick(rasterfile)
  Raster_QC   <- rast(rasterfile)[[12]]
  Raster_good <- mask(rasterfile, mask = Raster_QC, inverse = TRUE, maskvalues = QC_good_value)
  names(Raster_good) <- band_names
  return(Raster_good)
}
