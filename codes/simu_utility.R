library(randomForest)
library(doFuture)
library(plyr)


transformBeta <- function(mu, var) {
  # beta: https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
  MA <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  MB <- MA * (1 / mu - 1)
  return(c(MA, MB))
}

transformGamma <- function(mu, var) {
  # beta: WIKI
  MB <- mu / var
  MA <- mu * MB
  return(c(MA, MB))
}

################################################################################################################
## Post-process
################################################################################################################

# oneSettingMultiSeeds is a list of seeds, each is a list: list("over", "list of thres")
# rs[["ours"]] is a list of ws

putOurs_into_oneSettingMultiSeeds <- function(oneSettingMultiSeeds,
                                              r, method =  "ours") {
  n.seed = length(oneSettingMultiSeeds)
  count = 1
  for (i in 1:n.seed) {
    r.ours.a.list.of.ws = list()
    for (j in 1:length(ws)) {
      r.ours.a.list.of.ws[[j]] = r[[count]]
      count = count + 1
    }
    oneSettingMultiSeeds[[i]][[method]] = r.ours.a.list.of.ws
  }
  return(oneSettingMultiSeeds)
}

#' prepareDataForSIR
#' transform the current simulated data to input for getSIR
#' 
prepareDataForSIR <- function(trajectories) {
  n_city = length(trajectories)
  t0 = length(trajectories[[1]]$action + 1)
  data_SIR = list()
  for (i in 1:n_city) {
    data_i = list()
    data_i[["action"]] = c(trajectories[[i]]$action, NA)
    data_i[["infected"]] = trajectories[[i]]$state[3, ]
    data_i[["date"]] = trajectories[[i]][["date"]]
    data_i[["gdp"]] = trajectories[[i]]$gdp
    data_i[["population"]] = trajectories[[i]]$popu
    data_SIR[[i]] = data_i
  }
  return(data_SIR)
}


#' plotParetoFrontierforSimu
#' plot the Pareto and points for each city, in the simulation experiments
#' @param RC_fixed a list of RC vectors

plotParetoFrontierforSimu <-function(ours, oracle, point_size = 1, line_size = 1 
                                     , add_sd
                                     , sd_ours
                                     , sd_thre_XN
                                     , sd_thre_T3
                                     , sd_thre_T2
                                     , r_thre_XN
                                     , r_thre_T3
                                     , r_thre_T2
                                     , RC_fixed, behav, city_name, origin_zero = F, inlog = F, echo = F){
  figure_size = point_size
  
  R = ours[1,]
  if (inlog)R = log10(R)
  C = ours[2,]
  citytype = c(rep(city_name, length(C)))
  dat = data.frame(R = R, C = C, citytype)
  theme_set(theme_bw())
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # cbp1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
  #           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  policies = c("Pareto", "mitigation", "suppresion", "behaviour", "threshold", "oracle")
  colors <- c( "Pareto" = cbp1[4], "oracle" = cbp1[8], "threshold" = cbp1[2], "RC_fixed"=cbp1[3],
              "mitigation" = cbp1[6], "suppresion" = cbp1[1], "behaviour" = cbp1[7])
  shapes <- c("Pareto" = 16, "oracle" = 1
              , "threshold" = 5,"mitigation"= 2, "suppresion" = 0
              , "behaviour" = 13)
  
  ##########################################################################
  # ours
  ##########################################################################
  # scale_color_manual(name="Action",labels=c(1,2,3), values=colors, drop=FALSE) +
  # scale_shape_manual(name="Action",labels=c(1,2,3), values=c(15, 16, 17), drop=FALSE) +
  # , color=Action, shape=Action
  dat$policy = "Pareto"
  gg = ggplot() +
    geom_point(data = dat, aes(x=R, y=C, colour=policy, shape=policy), size = figure_size * 1.3) + # 1.3
    geom_line(data = dat, aes(x=R, y=C, colour=policy, shape=policy), show.legend = FALSE, alpha = 1, size = line_size) +
    xlab("Count of cumulative infected (in log scale)") + ylab("Loss of GDP (billion Chinese Yuan)")+
    # ggtitle(city_name) +
    facet_wrap(~citytype, scale="free_y", dir="v") +
    scale_colour_manual(name="Policy",values = colors, drop=F)+
    scale_shape_manual(name="Policy",values = shapes, drop=F)+
    scale_fill_manual(name="Policy",values = colors, drop=F)+
    labs(colour  = "Policy", shape = "Policy", fill="Policy", drop=F) +
    # theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.caption = element_text(color = "#D6604D", face = "italic")   ) 
  if(add_sd){gg = gg + geom_rect(data=dat, aes(xmin=R-sd_ours[1,]*1.96, xmax=R+sd_ours[1,]*1.96, ymin=C-sd_ours[2,]*1.96, 
                                     ymax=C+sd_ours[2,]*1.96, fill=policy, colour=policy),  alpha=0.1) }
  if(echo){print(gg)}

  ##########################################################################
  # competing
  ##########################################################################
 
  if (!missing(r_thre_T2)) {
    R_thre_T2 = r_thre_T2[1,]
    if (inlog)R_thre_T2 = log10(R_thre_T2)
    C_thre_T2 = r_thre_T2[2,]
    dat_thre_T2 = data.frame(R = R_thre_T2, C = C_thre_T2)
    dat_thre_T2$policy = "mitigation"
    dat_thre_T2$sdk_x = sd_thre_T2[1,]*1.96;
    dat_thre_T2$sdk_y = sd_thre_T2[2,]*1.96;
    gg =   gg + geom_point(data = dat_thre_T2, aes(x = R, y = C, color= policy, shape=policy), size=figure_size) + 
      geom_line(data = dat_thre_T2, aes(x = R, y = C, color= policy, shape=policy)
                , alpha = 1, show.legend = FALSE
                , size = line_size) 
    if(add_sd){
      gg = gg + geom_rect(data = dat_thre_T2, aes(xmin=R-sdk_x, xmax=R+sdk_x, ymin=C-sdk_y, 
                                                      ymax=C+sdk_y, fill=policy, colour=policy),  alpha=0.1)
    }
    if(echo){print(gg)}
  }
  if (!missing(r_thre_T3)) {
    R_thre_T3 = r_thre_T3[1,]
    if (inlog)R_thre_T3 = log10(R_thre_T3)
    C_thre_T3 = r_thre_T3[2,]
    dat_thre_T3 = data.frame(R = R_thre_T3, C = C_thre_T3)
    dat_thre_T3$policy = "suppresion"
    gg = gg + 
      geom_point(data = dat_thre_T3, aes(x = R, y = C, color= policy, shape=policy), size=figure_size) + 
      geom_line(data = dat_thre_T3, aes(x = R, y = C, color= policy, shape=policy)
                , alpha = 1, show.legend = FALSE
                , size = line_size)
    if(add_sd){
      gg = gg + geom_rect(data = dat_thre_T3, aes(xmin=R-sd_thre_T3[1,]*1.96, xmax=R+sd_thre_T3[1,]*1.96, ymin=C-sd_thre_T3[2,]*1.96, 
                                                        ymax=C+sd_thre_T3[2,]*1.96, fill=policy, colour=policy),  alpha=0.1)
    }
    if(echo){print(gg)}
  }
  
  if (!missing(behav)) {
    if(inlog) behav[1] = -log10(-behav[1])
    dat_behav = data.frame(R=-behav[1], C=behav[2], policy="behaviour")
    if(add_sd){
      gg = gg + geom_point(dat_behav, mapping = aes( x = R, y = C , color= policy, shape=policy, fill=policy), size= figure_size + 1) # , shape = 5
    } else{
      gg = gg + geom_point(dat_behav, mapping = aes( x = R, y = C , color= policy, shape=policy), size= figure_size + 1) # , shape = 5
    }
    if(echo){print(gg)}
  }  
  # dat_behav = data.frame(x =  -behav[1], y = behav[2])
  # gg = gg + geom_point(data = dat_behav, aes( x , y , color= policies[4], shape=policies[4]), size=figure_size) # , shape = 5
  
  if (!missing(r_thre_XN)) {
    R_thre_XN = r_thre_XN[1,]
    if (inlog)R_thre_XN = log10(R_thre_XN)
    C_thre_XN = r_thre_XN[2,]
    dat_thre_XN= data.frame(R = R_thre_XN, C = C_thre_XN)
    dat_thre_XN$policy = "threshold"
    gg = gg + 
      geom_point(data = dat_thre_XN, aes(x = R, y = C
                                         , color = policy, shape = policy) , size=figure_size) + 
      geom_line(data = dat_thre_XN, aes(x = R, y = C, color = policy, shape = policy)
                ,  alpha = 1, show.legend = FALSE
                , size = line_size)
    if (add_sd){
      gg = gg + geom_rect(data = dat_thre_XN, aes(xmin=R-sd_thre_XN[1,]*1.96, xmax=R+sd_thre_XN[1,]*1.96, ymin=C-sd_thre_XN[2,]*1.96, 
                                                      ymax=C+sd_thre_XN[2,]*1.96, fill=policy, colour=policy),  alpha=0.1)
    }
    
    if(echo){print(gg)}
  }
  

  ##########################################################################
  # limit
  ##########################################################################
  if(echo){print(gg)}
  xs = c(R, -behav[1], R_thre_XN, R_thre_T3, R_thre_T2)
  ys = c(C, behav[2], C_thre_XN, C_thre_T3, C_thre_T2)
  sds = cbind(sd_ours, c(0,0), sd_thre_XN, sd_thre_T3, sd_thre_T2) *2
  sdsx = sds[1,];  sdsy = sds[2,]
  const = 1
  xsmin = const * ( min(xs) - max(sdsx) )
  xsmax = const * ( max(xs) + max(sdsx) )
  ysmin = const * (min(ys) - max(sdsy))
  ysmax = const * (max(ys) + max(sdsy))
  if(origin_zero){
    gg = gg + xlim(0, xsmax) +
      ylim(0, ysmax)
    
  }else{
    gg= gg + xlim(xsmin, xsmax) +  ylim(ysmin, ysmax)
  }
  gg = gg +  theme(axis.title=element_blank()) + expand_limits(x=c(1,4))
  return(gg)
}



########################################################################################################
########################################################################################################

#' reportMeanStd
#' 
#' @return a len-rep list of RC = matrix(NA, 2, n_city)
#' @return list(mean, std) of RC = matrix(NA, 2, n_city)
#' 
reportMeanStd <-function(a, inlog = T){
  res.reps = list()
  for (i in 1:length(a)) {
    res.reps[[i]] = a[[i]][["RC"]]
  }
  n_city = dim(res.reps[[1]])[2]
  rep = length(res.reps)
  u_RC = Reduce('+', res.reps) / length(res.reps)
  sd_RC = matrix(NA, 2, n_city)
  
  for (i in 1:2) {
    for (j in 1:n_city) {
      value = rep(NA, rep)
      for (k in 1:rep) {
        value[k]  = res.reps[[k]][i,j]
      }
      if (inlog) {
        if (i == 1){
          sd_RC[i, j] = sd(log10(value))
        } else{
          sd_RC[i, j] = sd((value))
        }
        
      }else{
        sd_RC[i, j] = sd((value))
      }
    }
  }
  return(list(u = u_RC, sd = sd_RC))
}


#' getParetoCities
#' input:  name -> a (len-rep) list of len-n_method lists. If Pareto, then a len_n_weight list of RC = matrix(NA, 2, n_city)
# 

getParetoCities <-function(method, oneSettingMultiSeeds, ws, sd_out=FALSE, inlog = T){
  r_pareto = list()
  rep = length(oneSettingMultiSeeds)
  for (j in 1:length(ws)) {
    r = list()
    for (i in 1:rep) {
      r[[i]] = oneSettingMultiSeeds[[i]][[method]][[j]]
    }
    r_pareto[[j]]  = reportMeanStd(r, inlog = inlog)
  }
  
  u_cities = list()
  for (l in 1:n_city) {
    u = matrix(NA, 2, length(ws))
    for (j in 1:length(ws)) {
      u[,j] = r_pareto[[j]]$u[, l]
    }
    u_cities[[l]] = u
  }
  
  sd_cities = list()
  for (l in 1:n_city) {
    sd = matrix(NA, 2, length(ws))
    for (j in 1:length(ws)) {
      sd[,j] = r_pareto[[j]]$sd[, l]
    }
    sd_cities[[l]] = sd
  }
  
  if(sd_out){
    return(list(u=u_cities, sd=sd_cities))
  }else{
    return(u_cities)
  }
  
}



getFixedPolicy <-function(method, oneSettingMultiSeeds){
  over = list()
  rep = length(oneSettingMultiSeeds)
  for (i in 1:rep) {
    over[[i]] = oneSettingMultiSeeds[[i]][[method]]
  }
  r_over = reportMeanStd(over)
  return(r_over)
}

evaluateBehavPolicySimu <-function(data, i, t0, T_end, lag){
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




#' regetModel
#' for some unknown reasons, the classification model can not be reloaded. retrained each time.
#' 
#' @return a classifier

# regetModel <-function(path = "/Users/mac/Google Drive/Confident/out/augmented_data.RDS"){
#   # model = readRDS("/Users/mac/Google Drive/Confident/out/classifier.RDS")
#   
#   dat = readRDS(path)
#   
#   library(lightgbm)
#   
#   set.seed(0)
#   
#   last_dim = dim(dat)[2]
#   n = dim(dat)[1]
#   test_index = sample(1:n, 5000)
#   
#   train <- as.matrix(dat[-test_index, ])
#   test <- as.matrix(dat[test_index, ])
#   dtrain <- lgb.Dataset(data = train[, 1L:(last_dim - 1)], label = train[, last_dim])
#   dtest <- lgb.Dataset.create.valid(dtrain, data = test[, 1L:(last_dim - 1)], label = test[, last_dim])
#   valids <- list(test = dtest)
#   
#   # Method 1 of training
#   params <- list(objective = "multiclass", metric = "multi_error", num_class = 3L)
#   model <- lgb.train(
#     params
#     , dtrain
#     , nrounds = 5000
#     , valids = valids
#     , max_depth = 25
#     , max_bin = 1000
#     , eval_freq = 50
#     , num_threads = 1
#     # , min_data = 5L
#     , learning_rate = 1e-2
#     , early_stopping_rounds = 500L
#   )
#   # testing results
#   my_preds <- predict(model, test[, 1L:(last_dim-1)], reshape = T)
#   my_pred = apply(my_preds,1,which.max) - 1
#   
#   test_label = test[, last_dim]
#   acc = sum(my_pred == test_label) / dim(test)[1]
#   cat(acc)
#   return(model)
# }



