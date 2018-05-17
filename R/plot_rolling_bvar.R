plot.rolling_bvar <- function(rolling_bvar){

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

  print(gg_ridge)
  print(gg_line)

  list(gg_ridge, gg_line)

}
