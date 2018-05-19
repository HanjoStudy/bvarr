#' Evaluate density forecast with point realisations
#'
#' @details crps_sample employs an empirical version of the quantile
#'   decomposition of the CRPS (Laio and Tamea, 2007) when using method = "edf".
#'   For method = "kde", it uses kernel density estimation using a Gaussian
#'   kernel. The logarithmic score always uses kernel density estimation. The
#'   bandwidth (bw) for kernel density estimation can be specified manually, in
#'   which case it must be a positive number. If bw == NULL, the bandwidth is
#'   selected using the core function bw.nrd. Numerical integration may speed up
#'   computation for crps_sample in case of large samples dat.
#'
#'   A lower score indicates a better forecast
#'
#' @param rolling_bvar object of class `rolling_bvar`, see \link[rolling_BVAR]{bvarr}
#' @param plot produce forecast evaluation plots
#'
#' @references
#'
#' Evaluating simulation based forecast distributions: Krueger, F., Lerch, S.,
#' Thorarinsdottir, T.L. and T. Gneiting (2016): 'Probabilistic forecasting and
#' comparative model assessment based on Markov Chain Monte Carlo output',
#' working paper, Heidelberg Institute for Theoretical Studies, available at
#' http://arxiv.org/abs/1608.06802.
#'
#' Empirical quantile decomposition of the CRPS: Laio, F. and S. Tamea (2007):
#' 'Verification tools for probabilistic forecasts of continuous hydrological
#' variables', Hydrology and Earth System Sciences, 11, 1267-1277.
#'
#' @examples
#' df <- ggplot2::economics[,c(1, 4:6)]
#' bvar_res <- rolling_BVAR(df, date, start = "2013-01-01", by = 1, fixedWindow = F,
#'                           n_draws = 2000,
#'                           lags = 4,
#'                           cores = 3)
#' evaluate_bvar(bvar_res)
#'
#' @return Value of the score
#' @export
#'
#'
evaluate_bvar <- function(rolling_bvar, plot = T){
  if(!"rolling_bvar" %in% class(rolling_bvar)) stop("Only rolling evaluation at moment. Object not of valid class rolling_bvar, see bvarr::rolling_BVAR")

  actuals <- rolling_bvar$actuals %>%
    gather(., variable, metric, -date) %>%
    mutate(date = as.character(date)) %>%
    filter(date %in% rolling_bvar$forecast_density$date) %>%
    group_by(variable, date) %>%
    mutate(grouped_id = row_number()) %>%
    spread(date, metric) %>%
    select(-grouped_id) %>%
    group_by(variable) %>%
    nest(.key = actual) %>%
    mutate(actual = map(actual, ~.x %>% t %>% c))

  df_mcmc <- rolling_bvar$forecast_density %>%
    gather(., variable, metric, -date) %>%
    group_by(variable, date) %>%
    mutate(grouped_id = row_number()) %>%
    spread(date, metric) %>%
    select(-grouped_id) %>%
    group_by(variable) %>%
    nest(.key = mcmc) %>%
    mutate(mcmc = map(mcmc, ~.x %>% t %>% as.matrix))

  forecast_dates <- unique(rolling_bvar$forecast_density$date)

  df_forecast <- df_mcmc %>%
    left_join(actuals, by = c("variable")) %>%
    mutate(dates = map(list(forecast_dates) %>%  rep(3), ~.x)) %>%
    mutate(CRPS = map2(mcmc, actual, ~ scoringRules::crps_sample(dat = .x, y = .y)),
           LS = map2(mcmc, actual, ~ scoringRules::logs_sample(dat = .x, y = .y)))

  df_forecast <- df_forecast %>%
    select(variable, dates, CRPS, LS) %>%
    unnest(.) %>%
    gather(., type, value, -c(variable, dates))

  if(plot){
    df_forecast %>%
      ggplot(., aes(x = dates, y = value, fill = type))+
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_wrap(~type + variable, scales = "free_x")

    df_forecast %>%
      ggplot(., aes(value, fill = type)) +
      geom_density() +
      facet_wrap(~type + variable, scales = "free_x")
  }

  return(df_forecast)
}
