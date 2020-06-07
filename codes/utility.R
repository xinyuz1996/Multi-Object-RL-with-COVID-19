## Packages ---------

# library(rstanarm)
library(abind)
library(crayon)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(matrixStats)
library(parallel)






# install.load::install_load("invgamma")

# Transition functions ----------------------------------------------------
# Y1 = g(X_lt1_old, A1, X_lt1, X_l, w)
g <-function(X_lt_old, A, X_lt, GDP, w, seed = 1){
  R = g_R(X_lt_old, X_lt)
  C = g_C(A, GDP, seed = seed)
  return( c(utility = R - w * C, Reward = R, Cost = C) )
}

#' @param X_lt (S, I, R)
#' negative
g_R <-function(X_lt_old, X_lt){
  return(X_lt[1] - X_lt_old[1])
}


# positive
g_C <-function(A, GDP = 3031.149 / 365, 
               vec_c = c(0, 0.443, 0.518), sigma_c = c(0, sqrt(0.0479), sqrt(0.0354)),
               deterministic = F, seed = 1,
               ratio = 0.1
               ){
  if (deterministic)return(vec_c[A] * GDP * ratio)
  set.seed(seed)
  return(rnorm(1, vec_c[A], sigma_c[A]) * GDP * ratio)
}

#' f
#' State transition function, sampled based on the SIR model with current parameter estimates
#'
#' @param theta (sigma_R, gamma, sigma_S, beta) or (gamma, beta)
#' @param A  A_{l,t} among {1,...,J}
#' @param XR_1 if specified, it indicts the XR at next time is observable and should be respected.
#' @param population population for region l
#' @export nextstate  c(XS=XS_new, XI=XI_new, XR=XR_new)
#' @export poisson  whether or not to use poisson

f <- function(X_lt, A, XR_1, theta, f1 = NA, X_l = NA, MAP.only = F, 
              seed = 1){
  
  set.seed(seed)
  # get current state variables and parameters
  J = 3
  XS = X_lt[1]
  XI = X_lt[2]
  XR = X_lt[3]
  popu = sum(X_lt)

  f1 = f2 = 1
  # no f1, f2 yet
  if (MAP.only) {
    gamma = theta[1]
    beta = theta[1 + A]
  }else{
    gamma = rbeta(1, theta[1], theta[2])
    beta = rgamma(1, theta[2 + A], theta[2 + J + A])
  }
  
    if (missing(XR_1)){
      XS_new = XS -  rpois(1, max(beta * f1 * XS * XI / popu, 0))
      e_R = rbinom(1, max(as.integer(XI), 0), gamma * f2)
      XR_new = XR + e_R

      XI_new = popu - XS_new - XR_new
    }else{# next XR is observable
      error = rpois(1, beta * f1 * XS * XI / popu)
      XS_new = XS -  error
      XR_new = XR_1
      XI_new = max(popu - XS_new - XR_new, 0) # the popu is not consistent in the data?
    }
  r = c(XS = as.integer(XS_new), XI = as.integer(XI_new), XR = as.integer(XR_new))
  names(r) = c("XS","XI","XR")
  return(r)
}


# Helper  funuctions ------------------------------------------------------


#' Smooth moving average for a series

SMA <- function(series, weight_center){
  l = length(series)
  new_series = c(series[1])
  w = (1 - weight_center) / 2
  for (i in 2:(l-1)) {
    new_series = c(new_series, w * series[i - 1] + w * series[i + 1] + weight_center * series[i])
  }
  new_series = c(new_series, series[l])
  new_series = as.integer(new_series)
  
  return(new_series)
}

#' shift_state
#'
#' @import data $date $infected (counts) $action
#' @param Xl n by p, the original data for l^th region.
#' @param lag =9 shifted time length.
#' @param MA if specified, then do moving average for state variables
#'                 and the weight of the middle point.
#'
#' @export Xl_shifted a data frame with (n+lag) by p. The last lag elements of XI and XS are NAs.
#' @export XS 1 - XI - XR
#' @export XI shifted infectious population ~= future eight days increasing infectious population
#' @export XR removed population ~= cumulative confirmed cases
#' @export A 1 to J.
#' Xl_s() = shift_state(Xl) # Xl:1/24-3/5 to Xl_s:1/16-2/26
shift_state <- function(Xl, popu, MA, lag = 9){
  ###
  # MA:
  dates = Xl$date
  beforedates = rev(seq(as.Date(dates[1]),  by="-1 day", length.out = lag+1))[1:lag]
  date = c(beforedates, dates)

  I_confirmed = Xl$infected
  if (!is.na(MA)) {
    I_confirmed = SMA(I_confirmed, MA)
  }
  I_confirmed = c(rep(0, lag), I_confirmed)
  XI = data.table::shift(I_confirmed, n = lag, fill = NA, type="lead") - I_confirmed
  XR = I_confirmed
  if (missing(popu)) {
    XS = 1 - XI - XR
  }else{
    XS = popu - XI - XR
  }

  #action_shift = data.table::shift(Xl$action, n = lag, fill = NA, type="lead")
  A = c(rep(1, lag), Xl$action)
  #if(times_population)return(list(XI, XR, XS, A))
  shift_data <- data.frame(date, XS, XI, XR, A)
  return(shift_data)
}


## newInfect: calculate new infected number as XR_t - XR_{t-1}
newInfect <- function(removed){
  r = removed - data.table::shift(removed, n = 1, fill = 0, type="lag")
  r[1] = r[2]
  return(r)
}


#################################### Helper Functions  ####################################

#' getValueforThisLambda
#' 
#' @param
#'
getValueforThisLambda <-function(vlam, XR_t0_lag, A_t0_1_lag, popu, 
                                 t0, T_end, theta, GDP_daily, w, M
                                 , lag = 9
                                 , init_seed = 1
                                 , decision.t = 7
                                 , MC = F
                                 , Y_smallest = -1e10){
  ##
  T_t0 = T_end - t0 + 1
  
  
  init_belief_state = generateBeliefState(theta, XR_t0_lag, c(A_t0_1_lag, rep(0,T_t0 - 1))[1:(1+(lag - 1))], popu)
  
  one_traj <- function(m){
    X_lt = init_belief_state
    init_seed_m  = init_seed * m
    states = array(NA, c(3, T_t0 + lag))
    states[3, 1:(lag + 1)] = XR_t0_lag
    As = c(A_t0_1_lag, rep(0,T_t0 - 1))
    states[, 1 + lag] = init_belief_state
    
    YRC_kk = rep(0, 3)
    V_kk = R_kk = C_kk = 0 # for each trajectory
    for (tt in 1: (T_t0 - 1)) {
      
      if(tt >= (lag + 1))if (states[2, tt] == 0 ) break # early stop to save computation
      
      
      if (tt %% decision.t == 1) { # decision point
          belife_state = generateBeliefState(theta, states[3, tt:(tt+lag)], As[tt:(tt+(lag - 1))], popu)
          A =  policy(belife_state, vlam)
      }else{
        A = A
      }
      
      As[tt + lag] = A
      X_lt_old = states[, tt + lag]
      X_lt = f(X_lt = X_lt_old, A = A, theta = theta, seed = init_seed_m)
      
      states[, tt + lag + 1] = X_lt
      YRC_one_step = g(X_lt_old, A, X_lt, GDP_daily, w, seed = init_seed_m) # , deterministic = T
      YRC_kk = YRC_kk + YRC_one_step
      init_seed_m = init_seed_m + 1
      
    }
    return(YRC_kk) #  / (T_t0 - 1) # more practical
  }
  
  every_epoch = 10
  epoch = M %/% every_epoch
  YRC_k_list = list()
  for (i in 1:epoch) { 
    YRC_k_list_i = lapply(c(((i - 1) * every_epoch): (i * every_epoch)), one_traj) 
    YRC_k_list_i = do.call(rbind, YRC_k_list_i)
    YRC_k_list[[i]] = YRC_k_list_i
    YRC_k = do.call(rbind, YRC_k_list)
    Y_mean = colMeans(YRC_k)[1]
    if (Y_mean < 2 * Y_smallest)break # for saving cost purpose
  }
  YRC_mean = colMeans(YRC_k)
  YRC_std = colSds(YRC_k) / sqrt(i * every_epoch)

  return(c(YRC_mean, YRC_std))
}

### generateBeliefState (according to Lihong Li)  -----
#' given {XR_t}_{t = t0-lag}^{t_0} and {A_t}_{t=t_0-lag}^{t_0-1}, generate XI_{t_0} and therefore get S_{t_0}, the input of SPSA. Based on S_{t_0}, we select A_{t_0}.
#' similar with function Prediction. based on f()
#'
#' @param theta (gamma, beta)
#' @param XRs
#' @param As A \in \{1, ..., J\}
#' @param popu
#'
#' @Output (X^S, X^I, X^R)


generateBeliefState <- function(theta, pre_XRs, pre_As, popu
                                , MAP.only = F, lag = 9){
  # gamma = theta[1]
  J = 3

    XIG = pre_XRs[(lag + 1)] - pre_XRs[1] # X^I_{generated}
    f1 = 1
    for (i in 1:lag) {
      if (MAP.only) {
        beta = theta[1 + pre_As[i]]
      }else{ # get expectation of the posterior
        MAJ = theta[2 + pre_As[i]]
        MBJ = theta[2 + J + pre_As[i]]
        beta = MAJ / MBJ
      }
      new_infected = f1 * beta * (popu - XIG - pre_XRs[i]) * XIG / popu
      # observed_removed = gamma * f2 * XIG 
      observed_removed = pre_XRs[i + 1] - pre_XRs[i]
      XIG = XIG + new_infected  - observed_removed # the expected (mode?)
      XIG = max(XIG, 0)
    }
    X_lt = c(popu - XIG - pre_XRs[(lag + 1)], XIG, pre_XRs[(lag + 1)])

  return(X_lt)
}


#################################### Helper Functions  ####################################
norm_vec <- function(x) sqrt(sum(x^2))

#' Select action according to lambda
#' [0, lam2), [lam2, lam3), [lam3, inf)
#' @param Xlt (X^S, X^I, X^R), where X^I is the number of latent infectious
#' @param vlam lam_2, ..., lam_J
#'
#' @Output: if XI < lam_k, return action k
policy <-function(Xlt, vlam){
  XI = Xlt[2]
  J_1 = length(vlam)
  for (i in 1:J_1) {
    if (XI < vlam[i]) {
      return(i)
    }
  }
  return(J_1 + 1)
}


### evaluateBehavPolicy
#' evaluate the R&C of the observed performance of local goverments.
#'
#'
#' o.w., X^R_{T + lag} != X^R_{T}. We assume the government will follow the recommended policy afterwards for fairness.
#'
#' @param
evaluateBehavPolicy <-function(data, i, t0, T_end, shift = T, lag = 9){
  data_i = data[[i]]
  popu = data_i$population[1]
  GDP_daily = (data_i$gdp / 1e9 / 365)[1]
  XR = data_i$infected
  As = data_i$action
  
  V_R = XR[t0 + lag] - XR[T_end]
  V_C = 0
  for (t in t0:T_end) {
    V_C = V_C + g_C(As[t], GDP_daily, deterministic = T)
  }
  return(c(V_R, V_C))
}
