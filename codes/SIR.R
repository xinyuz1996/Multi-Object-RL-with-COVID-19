library(ggthemes)


#' f1_func
#' Heterogeneity function control transmission rate \beta
#'
#' @param data data as a list with L data frames for each city
#' @param location integer
f1_func <- function(data=data, location=1){
  return(1)
}

#' f2_func
#' Heterogeneity function control recover rate \gamma
#'
#' @param data data as a list with L data frames for each city
#' @param location integer
f2_func <- function(data=data, location=1){
  return(1)
}


#' getSIRMAP
#' index: [1,lag] + [lag + 1, lag + len(data))],  where (lag + 1):(lag + len(data))) represent data on 1:120.
#' Use observations until A_{t0-1 + 8} and X_{t0 + 8} to train the parameter theta.
#' @param data data as a list with L data frames for each city, data$infected, $suspected, $removed, $action, $population
#' @param t0 t0 <= dim(data)[0]. Use state until A_{t0-1} and S_{t0} to train the parameters, i.e., data until t0 + 8
#' @param lag shifted time between data and self defined state variable
#' @param J integer, number of action levels
#' @param u = c(u_1, ..., u_J) A vector of length J, j^th element representing prior for relative transmission rate of action j to action 1
#' @param region_index vector of integers, indicating which region you are gonna select for estimation
#' @param alpha the significance level used to generate prediction interval for estimated parameters
#' @param MA = 0.4 central weight for smooth moving average method
#' @param preprocess_tl whether to use redfined tl, instead of the first non-zero
#' @param echo boolean variable indicating if the zlts/zltr should be output during the process
getSIRMAP <- function(data=data, t0=53, lag=9,J=3, u=c(1, .2, .1), region_index=c(1:6), alpha=0.1, 
                      MA_index = NA, preprocess_tl = T, echo = F, brief = T, MAP.only = F, 
                      staticgamma = T, for_simu = F){

  if (for_simu) {
    t_max = length(data[[region_index[1]]]$infected)
  }else{
    t_max = dim(data[[region_index[1]]])[1]
  }
  
  if((t0) > t_max){
    stop(paste("The input t0 should be an integer from 1 to", t_max,".\n"))
  }
  
  XRSUM <- c(0)
  XISUM <- c(0)
  XSSUM <- rep(0,J)
  XSISUM <- rep(0,J) # XS*XI/Ml
  ASUM <- rep(0,J) # counts of action
  
  ## collect transition data  
  for (l in region_index){
    f1 = f1_func(data, l)
    # f2 = f2_func(data, l)
    Xl = addTltoData(data[[l]], preprocess = preprocess_tl)
    tl = Xl$tl[1]
    Ml = Xl$population[1] # supposed to be variant, but invariant based on data collection
    # If there's no available data, then skip this city
    if(is.na(tl)){
      next
    }
    
    # 9:50 shift to 1:42
    Xl_shifted = shift_state(Xl, popu=Ml, lag = lag, MA = F)
    XS = as.integer( Xl_shifted$XS )
    XI = as.integer( Xl_shifted$XI )
    XR = as.integer( Xl_shifted$XR )
    A =  Xl_shifted$A
    
    if((t0-1) >= (max(tl, 0))){
      for (t in (max(tl, 0)):(t0-1)){
        
        Xti = XI[t]
        Xtr = XR[t]
        Xts = XS[t]
        
        Xti_p1 = XI[t+1]
        Xtr_p1 = XR[t+1]
        Xts_p1 = XS[t+1]
        
        if ( Xti > 0 ){
          
          # update for invariant gamma
          Zltr =  (Xtr_p1 - Xtr)
          XRSUM = XRSUM + Zltr
          XISUM = XISUM + Xti
          
          # Update for action-varying beta
          Zlts =  (Xts - Xts_p1 )
          Zltsi = (Xts/Ml)*Xti*f1 # Numerical stable
          
          i = A[t]
          ASUM[i] = ASUM[i] + 1 # counter
          XSSUM[i] = XSSUM[i] + Zlts
          XSISUM[i] = XSISUM[i] + Zltsi
          #if (i == 1) {
          #  cat(XSSUM[i],Zlts,Xts,Xts_p1, "!!\n")
          #}
          
          if(echo){cat("l:",l,"t",t,"xti",Xti,": Zltr=",Zltr,"  .","Zlts=", Zlts,  "Zltsi=", Zltsi, "action=", i,"\n")}
        }
      }
    }
    
  }
  
  if(echo){cat("ASUM:",ASUM)}
  
  
  ### Priors
  
  # Manipulate priors based on SARS
  R0_1.prior = 3.1511
  gamma.prior = 0.0821 
  beta_1.prior = R0_1.prior * gamma.prior 
  var.gamma.prior = .1
  
  # gamma
  bR = 2 * 1e3 # 1e3 average sum XI
  aR = gamma.prior/(1 - gamma.prior) * bR

  MA = aR + XRSUM
  MB = bR - XRSUM + XISUM
  
  gamma.mean = (MA) / (MB + MA)
  gamma.mode = (MA-1) / (MA+MB-2)
  gamma.var = (MA*MB)/( (MA+MB)^2*(MA+MB+1) )
  gamma.std = sqrt(gamma.var)
  gamma.95cilb = qbeta(c(alpha/2), shape1 = MA, shape2 = MB)
  gamma.95ciub = qbeta(c(1-alpha/2), shape1 = MA, shape2 = MB)# * sqrt(gamma.sigma2) + gamma.mu
  gamma = cbind(mean = gamma.mean, std = gamma.std, CILB = gamma.95cilb, CIUB = gamma.95ciub)
  row.names(gamma) = "gamma"
  
  
  ## betas
  bS = rep(2, J) * 1e3
  # beta1
  beta1.mean = (beta_1.prior * bS[1] + XSSUM[1]) / (bS[1] + XSISUM[1])
  
  beta.prior = c(beta_1.prior, beta_1.prior * u[2:J])

  aS = beta.prior * bS
  # other beta's
  
  
  
  
  # Posterior mean
  MAJ = aS + XSSUM
  MBJ = bS + XSISUM
  if(echo){cat(sep=" ","XRSUM:",XRSUM, ";XISUM=",XISUM, ";XSSUM",XSSUM,";XSISUM",XSISUM,"\n.")}
  
  
  
  # beta
  beta.mean =  MAJ / MBJ
  beta.mode = (MAJ-1) / MBJ
  beta.var = MAJ / (MBJ)^2
  beta.std = sqrt(beta.var)
  beta.95cilb = qgamma(c(alpha/2), shape = MAJ, rate = MBJ)
  beta.95ciub = qgamma(c(1-alpha/2), shape = MAJ, rate = MBJ)
  beta = cbind(mean = beta.mean, std = beta.std, CILB = beta.95cilb, CIUB = beta.95ciub)
  row.names(beta)=paste("beta",c(1:J), sep="")
  
  # R0 new
  R0 = matrix(nrow=3, ncol=4)
  nrep=1e5
  betas = matrix(rgamma(3*nrep, shape = MAJ, rate = MBJ), ncol=3, byrow = T)
  gammas = matrix( rbeta(nrep, shape1 = MA, shape2 = MB), ncol=1, byrow = T )
  for (i in 1:3){
    R0i = betas[,i]/gammas
    R0[i,] = c(mean(R0i), sd(R0i), quantile(R0i, c(.05, .95)))
  }
  colnames(R0) <- c("mean","std","CILB","CIUB")
  row.names(R0) = paste("R0",c(1:J), sep="")
  
  
  # combine parameters
  parameters = as.data.frame(round( rbind(gamma, beta, R0), 4 )) 
  sirMAP = round(parameters$mean, 3)
  sir.posterior = c(MA, MB, MAJ, MBJ, sirMAP)
  
  if (brief){
    if (MAP.only) {
      return(sirMAP)
    }else{
      return(sir.posterior)
    }    
  } else {
    return(parameters)
  }
}




#' addTltoData
#' add the tl to input data
#'
#' @param city_data data frames for each city, which is data[[l]]
addTltoData <- function(city_data, preprocess = FALSE){
  CI <-  city_data$infected
  n <- length(CI)
  i=0
  if (preprocess) {#' Two days with new confirmed in a row
    while ( (i <= (n-3))){
      if( (CI[i+1] < CI[i+2]) && (CI[i+2] < CI[i+3]))break
      i=i+1
    }
  }
  # consider case when there's no input WV us[[37]]
  while ( (i <= (n-1) && (CI[i+1]==0)  )){
    i=i+1
  }
  if (i == n){
    #cat("This city is safe by now, without any cases confirmed.")
    tl = NA
  } else{
    i = i+1
    tl = rep(i, n)
  }
  city_data$tl = tl
  return(city_data)
}


#' Prediction
#' Make prediction based on given theta
#' index: [1,lag] + [lag + 1, lag + len(data))]
#'
#' @param data  data_cn or data_us
#' @param theta parameter Prediction
#' @param l index of city between 1:length(data)
#' @param t_begin = 21 the begin date for prediction, 0 < t_begin <= data_end
#' @param t_end = 53 the end day for prediction date, t_begin<=t_end, no upper limit.
#' @param rep=1000 number of replicates to get the MC PI for future predictions. The larger, the smoother.
#' @param alpha = 0.05, significance level = 100(1-alpha)%
#' @param actions if specified, then it is a vector of action levels at {1, 2, ..., J} with length = t_end-t_begin+1 corresponding to t_begin:t_end with index (+8)
#' @param MA = NA central weight for smooth moving average method, between 0 to 1. Can not be used with poisson
#'
#' @export newstate = list(XSest=XSest * population, XIest = XIest * population, XRest = XRest * population, XNest = XNest * population)
#'
prediction <- function(data, theta, l, actions=NA, t_begin = 26, t_end = 53,
                       rep=1e3, alpha = 0.1,  MA_index = NA, echo=FALSE, lag=9, 
                       seed.method = "more", MAP.only = F
){
  ## exceptional cases 
  L = length(data)
  if((l<0) | (l>L)){
    stop(paste("The input city index l should be an integer from 1 to", L,".\n"))
  }
  if((t_begin<0)){
    stop(paste("Please input a positive t_begin larger than 0.\n"))
  }
  if((t_end< t_begin)){
    stop(paste("Please input a t_end no less than t_begin.\n"))
  }

  #' index: [1,lag] + [lag + 1, lag + len(data))]
  preprocess_tl = T
  
  ## preprocess data 
  f1 = f1_func(data, l)
  f2 = f2_func(data, l)
  Xl = addTltoData(data[[l]], preprocess = preprocess_tl)
  population = Xl$population[1]
  
  # shifted states
  Xl_shifted = shift_state(Xl, MA = MA_index, popu = population, lag=lag)
  XS = Xl_shifted$XS
  XI = Xl_shifted$XI
  XR = Xl_shifted$XR
  A =  Xl_shifted$A
  X_lts = cbind(XS, XI, XR) # why keep the part before?
  
  ## Initialization 
  XIest = matrix(0,nrow = c(t_end-t_begin+1), ncol=3)
  rownames(XIest) = c(t_begin:t_end)
  XNest = XSest = XRest = XIest
  XNtraj = XStraj = XRtraj = XItraj = matrix(0, nrow = rep, ncol=c(t_end-t_begin+1))
  
  ## 
  if (t_end > dim(Xl)[1] + lag + 1){ 
    X_lts = rbind(X_lts, matrix(0, nrow=(t_end-dim(Xl)[1]+1), ncol=3  ))
  }else{
    X_lts = rbind(X_lts, rep(0,3))
  }

  # Generate trajectories -----
  for (r in 1:rep){
    count = 1
    for (t in t_begin:t_end){
      
      if (seed.method == "more") {
        seed = r * 123 + count
      }else{
        seed = r  
      }
      
      if(count <= lag){ # XR is observable
        if(is.na(actions)){
          X_lts[t,] = f(X_lt = X_lts[t - 1,], A = A[t - 1], theta = theta, XR_1 = X_lts[t, 3], seed = seed, MAP.only = MAP.only)
          if (echo){cat("X_lts[t,]",X_lts[t,])}
        }else{
          X_lts[t,] = f(X_lts[t-1,], actions[t - t_begin + 1],theta = theta, XR_1 = X_lts[t, 3], seed = seed, MAP.only = MAP.only)
        }
      }else{ # XR is NOT observable
        if(is.na(actions)){
          X_lts[t,] = f(X_lts[t-1,], A[t-1], theta = theta, seed = seed, MAP.only = MAP.only)
          if (echo){cat("X_lts[t,]",X_lts[t,],"\n")}
        }else{
          X_lts[t,] = f(X_lts[t-1,], actions[t - t_begin + 1], theta=theta, seed = seed, MAP.only = MAP.only)
        }
      }
      count = count + 1
      
    }
    
    if (is.na(sum(X_lts[t,]))) {
      stop("is.na(sum(X_lts[t,]))")
    }
    
    XStraj[r,] = X_lts[t_begin:t_end,1]
    XItraj[r,] = X_lts[t_begin:t_end,2]
    XRtraj[r,] = X_lts[t_begin:t_end,3]
    XNtraj[r,] = newInfect(X_lts[t_begin:t_end,3])
  }
  
  ## Prepare output: calculate means and PIs ----
  XSest = t( apply(XStraj, 2, quantile, c(alpha / 2, .5, 1 - alpha / 2)) )
  XSest[,2] = colMeans(XStraj); colnames(XSest)[2] = "mean"
  XIest = t( apply(XItraj, 2, quantile, c(alpha / 2, .5, 1 - alpha / 2)) )
  XIest[,2] = colMeans(XItraj); colnames(XIest)[2] = "mean"
  XRest = t( apply(XRtraj, 2, quantile, c(alpha / 2, .5, 1 - alpha / 2)) )
  XRest[,2] = colMeans(XRtraj); colnames(XRest)[2] = "mean"
  XNest = t( apply(XNtraj, 2, quantile, c(alpha / 2, .5, 1 - alpha / 2)) )
  XNest[,2] = colMeans(XNtraj); colnames(XNest)[2] = "mean"
  
  rownames(XSest) = c(t_begin:t_end)
  rownames(XIest) = c(t_begin:t_end)
  rownames(XRest) = c(t_begin:t_end)
  rownames(XNest) = c(t_begin:t_end)
  newstate = list(XSest=XSest , XIest = XIest , XRest = XRest , XNest = XNest )
  
  return(newstate)
}


#' plotNewValidation
#' Validation for chinese data
#' index: [1,8] + [9, 8 + len(data))]
#' Trained by data from 1 to t0.
#' Test on data from t0 + 9 to t_end = 8 + len(data)
plotNewValidation <- function(data=data_cn, l, t0=t0, t_end = 61, alpha=.05, theta=theta, lag=9, 
                              MA_index = NA, MAP.only = F, onlyXR=T,
                              seed.method ="more", date_break="2 week"){

  
  names = names(data)
  ## Do preditction
  t_begin = t0 + 1
  newstate = prediction(data=data, theta=theta, l=l, t_begin = t_begin
                        , seed.method = seed.method
                        , t_end=t_end, alpha = alpha, MA_index = MA_index, lag = lag, MAP.only = MAP.only)
  
  ## Prepare observed data
  Xl = data[[l]]
  T_data = dim(Xl)[1]
  if(t_end > T_data+lag){
    stop(paste("The input t_end should be an integer from 1 to", T_data+lag,".\n"))
  }
  
  name = names(data)
  population = Xl$population[1]
  f1 = f1_func(data, l)
  f2 = f2_func(data, l)
  
  Xl_shifted = shift_state(Xl, MA = MA_index, popu = population, lag = lag)
  XI = Xl_shifted$XI
  XR = Xl_shifted$XR
  XS = Xl_shifted$XS
  
  ## Date indexes and record range
  xlab = as.Date(as.character(Xl$date), "%Y-%m-%d")
  x_lab_all = c(as.Date(as.character(Xl$date), "%Y-%m-%d")[1:lag] - lag, as.Date(as.character(Xl$date), "%Y-%m-%d"))
  x_lab_shifted = as.Date(as.character(Xl$date), "%Y-%m-%d") - lag
  
  x_range = c(1:(t_end))
  x_range_all = c(1:(t_end))
  # + 1 due to the as.POSIXct below
  x_date = as.POSIXct(xlab + 1)
  x_date_all = as.POSIXct(x_lab_all + 1)
  x_date_shifted = as.POSIXct(x_lab_shifted + 1)
  
  
  # align date for Observation state variable
  XR = data.frame(x_date = x_date_all[1:t_end], count = XR[1:t_end])
  XNew = data.frame(x_date = x_date_all[1:t_end], count = newInfect(XR[1:t_end,2]))
  # combine the starting point for Prediction state variable
  XRest = rbind(rep(XR[t_begin-1+lag, 2], 3), newstate$XRest[(1+lag):dim(newstate$XRest)[1],])
  XNest = rbind(rep(XNew[t_begin-1+lag, 2], 3), newstate$XNest[(1+lag):dim(newstate$XNest)[1],])
  
  # dimension
  main_title_suffix = ""
  colours = c("black","red", "red","orangered", rgb(0.2,0.4,0.1,0.7), rgb(0.8,0.4,0.1,0.7))
  main_title = paste("Predict by first", as.character(t_begin - 1), "days data", main_title_suffix)
  
  # date index for Prediction
  x_pred_date_R = x_date_all[(t_begin-1+lag):t_end]
  XRest = data.frame(x_date = x_pred_date_R, count = XRest)
  XNest = data.frame(x_date = x_pred_date_R, count = XNest)
  names(XRest) = c("x_date", "lower", "mean", "upper")
  names(XNest) = c("x_date", "lower", "mean", "upper")
  
  
  title_prefix=""
  
  
  if(onlyXR){
    type = c(rep(1,dim(XR)[1] ) )
    type = factor(type, levels=c(1,3,4), labels=c("Total confirmed cases","95%CI","4"))
    citytype = c(rep(names[l], dim(XR)[1]))
    dataobs = cbind(rbind(XR), type, citytype)
    names(dataobs) = c("date", "count","type","citytype")
    
    esttype = c(rep(1,dim(XRest)[1] ))
    citytype = c(rep(names[l], dim(XRest)[1]))
    esttype = factor(esttype, levels=c(1,3,4), labels=c("Total confirmed cases","95%CI","4"))
    dataest = cbind(rbind(XRest), type=esttype, citytype)
    colors <- c("Observation" = "#808080", "95%Bounds" = "orange", "Prediction" = "red","Action"="#0066CC")
    lines <- c("Observation" = 1, "95%Bounds" = 3, "Prediction" = 5, "Action"=4)
    gg = ggplot(dataobs) +
      theme_bw() +
      facet_wrap(~citytype, scale="free_y", dir="v") +
      geom_line(aes(date, count, color="Observation", linetype="Observation"), na.rm = T) +
      geom_line(data = dataest, aes(x = x_date, y = lower, color="95%Bounds",linetype = '95%Bounds')) +
      geom_line(data = dataest,aes(x = x_date, y = upper, color="95%Bounds", linetype = '95%Bounds')) +
      geom_line(data = dataest,aes(x = x_date, y = mean, color="Prediction",linetype="Prediction")) +
      geom_vline(xintercept = x_date_all[t_begin - 1 + lag],color="grey", linetype="longdash", size=.1) +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_datetime(date_breaks = date_break, date_labels = ("%b %d"))+
      labs(x="Date",y="count",color  = "Lines", linetype = "Lines", shape = "Lines", drop=F) +
      scale_color_manual(values = colors, drop=F) +
      scale_linetype_manual(values = lines, drop=F) +
      geom_ribbon(data=dataest, aes(x = x_date, ymin=lower,ymax=upper), fill=colors['95%Bounds'],alpha=.2)
  }
  date_range = (range( dataest[2:(t_end-lag-t0+1),1] ))
  deviations = dataest[2:(t_end-lag-t0+1),'mean'] - dataobs[(t0+lag+1):t_end,'count']
  RMSE = sqrt( sum(deviations^2 / length(deviations)) )
  return(list(gg = gg, 
              dataobs = dataobs, 
              dataest = dataest,
              deviations = deviations,
              RMSE = RMSE, 
              date_range = date_range))
}
