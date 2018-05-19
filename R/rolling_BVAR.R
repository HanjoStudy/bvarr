#' @title Rolling/expanding Bayesian VAR estimate
#'
#' @description Performs a rolling/expanding estimation of BVAR forecast estimates with re-estimation at each point t. Later implementation will use a model opbject
#' @param df data.frame containg Y variables and a date column
#' @param date_col specify column name which contains t time elements
#' @param start start date from which the rolling estimation should begin ex. "2010-01-01"
#' @param h what should the forecast horison be (Now only h = 1, by = 1)
#' @param fixedWindow whether it should be expaning or rolling. Fixed = T is rolling
#' @param include keep raw for sample analysis
#' @param n_draws number of draws to keep
#' @param lags number of lags for VAR estimate
#' @param cores if cores are specified then forecasts are implemented in parallel
#' @examples
#'
#' df <- ggplot2::economics[,c(1, 4:6)]
#'
#'
#' bvar_res <- rolling_BVAR(df, date, start = "2013-01-01", by = 1, fixedWindow = F,
#'                   n_draws = 2000,
#'                   lags = 4,
#'                   cores = 3)
#'
#' bvar_res$forecast_density %>% filter(date == "2015-02-01") %>%
#'   tidyr::gather(., metric, value, -date) %>%
#'   ggplot(., aes(value)) +
#'   geom_density() +
#'   facet_wrap(~metric, scales = c("free"))
#'
#' plot(bvar_res)
#'
#' @return data.frame
#' @export
#'

rolling_BVAR <- function(df, date_col, start, h = 1, fixedWindow = T,
                         include = "raw",
                         n_draws = 2000,
                         lags = 4,
                         cores = 2){

  options(stringsAsFactors = F)
  if(missing(date_col)) stop("No appropriate date column found: please provide data.frame with aptly named column: Date or date")
  if(missing(cores)) parallel <- FALSE
  # Prepare data slices
  orig <- df %>%
    mutate_if(is.factor, as.character)
  date_col <- dplyr::enquo(date_col)

  dates <- df %>%
    dplyr::pull(!!date_col)

  df <- df %>%
    dplyr::select(-!!date_col)

  initialWindow <- which(dates == start) - 1

  create_Time_slice <- function (y, initialWindow, horizon = 1, fixedWindow = fixedWindow,
            skip = 0)
  {
    stops <- seq(initialWindow, (length(y) - horizon), by = skip + 1)
    if (fixedWindow) {
      starts <- stops - initialWindow + 1
    }
    else {
      starts <- rep(1, length(stops))
    }
    train <- mapply(seq, starts, stops, SIMPLIFY = FALSE)
    nums <- gsub(" ", "0", format(stops))
    names(train) <- paste("Training", nums, sep = "")
    train
  }

  TS_slices <- create_Time_slice(y = 1:nrow(df), initialWindow, horizon = 1, fixedWindow = F, skip = 0)

  names(TS_slices) <- dates[-c(1:initialWindow)]

  # Run BVAR
  include <- ifelse(include == "raw", "raw", c("mean", "median", "sd", "interval", "raw"))

  if(parallel == T)
  {

    library(doSNOW)

    cl <- makeSOCKcluster(cores)

    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
    clusterExport(cl = cl, c('TS_slices', 'df'), envir = environment())

    pb <- txtProgressBar(min = 1, max = length(TS_slices), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)

    opts <- list(progress=progress)

    BVAR_forecast <- foreach(i = 1:length(TS_slices),
                       .packages = c("bvarr", "dplyr"),
                       .options.snow = opts) %dopar% {
                         x <- df[TS_slices[[i]],] %>% data.frame
                         setup <- bvar_conj_setup(x, p = lags)
                         model <- bvar_conj_estimate(setup, n_draws, verbose = T)
                         y <- bvar_conj_forecast(model, Y_in = NULL, Z_f = NULL, output = c("long", "wide"), h = h,
                                                                  out_of_sample = TRUE, type = c("prediction", "credible"), level = c(80, 95),
                                                                  include = include, fast_forecast = FALSE,
                                                                  verbose = FALSE) %>%
                           data.frame(date = names(TS_slices)[i], .)
                         if(include == "raw"){
                           if(h == 1){
                             y <- y %>%
                               purrr::set_names(c("date", names(x)))
                           }else if(h > 1){
                             y <- y %>%
                               data.frame %>%
                               select(1, tail(names(.), ncol(df))) %>%
                               purrr::set_names(c("date", paste0(names(df), paste0("(t+", h ,")"))))
                           }
                         }
                         y
                       }
    close(pb)

    BVAR_forecast <- do.call(rbind, BVAR_forecast)

  } else {

    BVAR_forecast <- list()
    for(i in 1:length(TS_slices)){
      x <- df[TS_slices[[i]],] %>% data.frame
      setup <- bvar_conj_setup(x, p = lags)
      model <- bvar_conj_estimate(setup, n_draws, verbose = T)
      y <- bvar_conj_forecast(model, Y_in = NULL, Z_f = NULL, output = c("long", "wide"), h = h,
        out_of_sample = TRUE, type = c("prediction", "credible"), level = c(80, 95),
        include = include, fast_forecast = FALSE,
        verbose = FALSE) %>%
        data.frame(date = names(TS_slices)[i], .)
      if(include == "raw"){
        if(h == 1){
          y <- y %>%
            purrr::set_names(c("date", names(x)))
        }else if(h > 1){
          y <- y %>%
            data.frame %>%
            select(1, tail(names(.), ncol(df))) %>%
            purrr::set_names(c("date", paste0(names(df), paste0("(t+", h ,")"))))
        }
      }
      BVAR_forecast[[i]] <- y
    }
    BVAR_forecast <- do.call(rbind, BVAR_forecast)
  }
  out <- BVAR_forecast %>%
    tbl_df %>%
    mutate_if(is.factor, as.character)

  out <- list(forecast_density = out, actuals = orig)
  class(out) <- c(class(out), "rolling_bvar")
  return(out)
}

#' @title Plot object of class rolling_bvar
#' @description plots the distribution and line plots for the forecasts using ggplot
#' @export
#'
#'
plot.rolling_bvar <- function(rolling_bvar, plot = T){

  if(!"rolling_bvar" %in% class(rolling_bvar)) stop("Object not of valid class rolling_bvar, see bvarr::rolling_BVAR")

  Yraw_test <- rolling_bvar$actuals %>%
    gather(., variable, metric, -date) %>%
    mutate(date = as.character(date)) %>%
    filter(date %in% rolling_bvar$forecast_density$date)

  gg_ridge <- rolling_bvar$forecast_density %>%
    gather(., variable, metric, -date) %>%
    ggplot(., aes(x = metric, y = date, group = date)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
    geom_point(aes(x = metric, date), color = "#fe0301" , size = 2, data = Yraw_test, inherit.aes = F) +
    facet_wrap(~variable, scales = "free_x") +
    theme_ridges() +
    labs(title = "Forecast density compared to actuals", subtitle ="Actuals (red)")


  gg_line <- rolling_bvar$forecast_density %>%
    gather(., variable, metric, -date) %>%
    group_by(variable, date) %>%
    summarise(median_est = median(metric)) %>%
    ungroup %>%
    left_join(Yraw_test, by = c("variable", "date")) %>%
    rename(., Estimates = median_est, Actuals = metric) %>%
    gather(., Type, Value, -c(variable, date)) %>%
    ggplot(., aes(date, Value, color = Type, group = Type)) +
    geom_line() +
    facet_wrap(~variable, scales = "free_y") +
    theme_ridges()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    labs(title = "Forecast point estimate compared to actuals", subtitle ="Actuals (red), Median Posterior density (blue)")
  if(plot){
    print(gg_ridge)
    print(gg_line)
  }

  return(list(gg_ridge, gg_line))
}


