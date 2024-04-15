library(zoo)

fit_loglinear = function(x, y, tries = 5){
  m1 <- NA
  #figure out start lag guess
  ymax = max(y)
  ymin = min(y)
  saturation_time = x[y > (ymin + ((ymax - ymin) * 0.95))][1]
  lag_time = x[y > (ymin + ((ymax - ymin) * 0.1))][1]
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ log_linear(x, r, y0, lag, end),
                        start = list(r = runif(1, 0.1, 2), 
                                     y0 = min(y), 
                                     lag = lag_time, 
                                     end = saturation_time)),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1)){
    return(c(NA, NA,NA,NA))
  }
  #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
  #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
  r = coef(m1)[1]
  y0 = coef(m1)[2]
  start = coef(m1)[3]
  end = coef(m1)[4]
  return(c(r, y0, start, end))
}


fit_logistic = function(x, y, tries = 5){
    m1 <- NA
    #figure out start lag guess
    lag_guess = guess_half_max(x,y)
    while(is.na(m1[1]) && (tries > 0)){
      m1 <- tryCatch( nls(y ~ logistic(x, r, lag, K, y0),
                          start = list(r = runif(1, 0.1, 1), 
                                       lag = lag_guess + sample(-20:20, 1), 
                                       K = max(y), 
                                       y0 = min(y))),
                      error = function(e) return(NA))
      tries = tries - 1
    }
    if (is.na(m1)){
      return(c(NA, NA,NA,NA))
    }
    #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
    #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
    r = coef(m1)[1]
    lag = coef(m1)[2]
    K = coef(m1)[3]
    y0 = coef(m1)[4]
    return(c(r, lag, K, y0))
}

fit_baranyi = function(x, y, tries = 100){
  m1 <- NA
  #figure out start lag guess
  lag_guess = guess_half_max(x,y)
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ baranyi(x, r, lag, ymax, y0),
                        start = list(r = runif(1, 0.1, 0.4), 
                                     lag = lag_guess + sample(-20:20, 1), 
                                     ymax = max(y), 
                                     y0 = min(y))),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1)){
    return(c(NA, NA,NA,NA))
  }
  #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
  #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
  r = coef(m1)[1]
  lag = coef(m1)[2]
  ymax = coef(m1)[3]
  y0 = coef(m1)[4]
  return(c(r, lag, ymax, y0))
}

guess_half_max = function(x, y){
  # this function looks at a logistically-increasing time series
  # and guesses when the growth rate is at a maximum (which is also when
  # the trend is at half-max)
  
  # put in order
  y = y[sort(x, index.return = TRUE)$ix]
  x = sort(x)
  
  #find approximate time of max growth rate, using diffs, on smoothed y
  y2 = diff(y)
  y2 = zoo::rollmean(y2, 11, fill = NA, align = "center")
  half_max_idx = which(y2 == max(y2, na.rm = TRUE))[1]
  half_max = x[half_max_idx]
  return(half_max)
}

baranyi <- function(t, r, lag, logymax, logy0){
  At = t + (1 / r) * log(exp(-r * t) + exp(-r * lag) - exp(-r * (t + lag)))
  logy = logy0 + r * At - log(1 + (exp(r * At) - 1) / exp(logymax - logy0))
  return(logy)
}

logistic <- function(t, r, lag, K, y0){
  y0 + (K - y0) / (1 + exp(-r * (t-lag)))
}

log_linear <- function(t, r, y0, lag, end){
  y = rep(0, length(t))
  y[t < lag] = y0
  y[t >= lag & t < end] = y0 + (t[t >= lag & t < end]-lag) * r
  y[t >= end] = y0 + (end - lag) * r
  return(y)
}
