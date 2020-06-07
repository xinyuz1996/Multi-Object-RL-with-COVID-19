
# Environment --------------------------------------------------
rm(list=ls())

## Packages
library(gtools)
library(dplyr)
library(tidyr)

## Paths
main_path = "./"
detectCores(); Sys.info()["sysname"]
data_path = paste(main_path,"data.rds",sep = "")
name = names(data)
code_path = paste(main_path,"code/code in package/",sep = "")
out_path = paste(main_path,"out/",sep = "")
figure_path = paste(main_path,"figure/final/",sep = "")
cv_path = paste(main_path,"out/cv/", sep="")

## Codes
miceadds::source.all(paste(code_path,sep = ""))

## Data
data = readRDS(data_path)
region_index = c(1:6)
name = names(data);name



######################################################################################################

# INitialize parameters
t0_train=61; t0_test = 12
alpha = .01; t_end = 120; lag=9


getRMSEforCitys <- function(train_region_index, test_region_index, t0_train = 61, t0_test = 12,
                            alpha = .1, t_end = 120, only.120 = T, ratio = T){
    theta = getSIRMAP(data=data, region_index = train_region_index, 
                    t0 = t0_train, )
  cat("R0:", theta[13:15],"\n")
  deviations <- c()
  for (l in test_region_index){
    p = plotNewValidation(data=data, l=l, t0=t0_test, theta=theta, t_end= t_end, alpha = alpha, lag = 9, date_break="30 day")
    if (only.120) {
      if (ratio) {
        deviations = c(deviations, p$deviations[length(p$deviations)] / p$dataobs$count[length(p$dataobs$count)])
      }else{
        deviations = c(deviations, p$deviations[length(p$deviations)])
      }
    }else{
      deviations = c(deviations, p$deviations)
    }
  }
  if (only.120) {
    return(deviations)
  }else{
    RMSE = sqrt( sum(deviations^2 / length(deviations)) )
    return(RMSE)
  }
  
}

# leave one ---------------------------------------------------------------

leaveOneCV <- function(i, only.120, ratio){
  train_region_index = c(1:6)[-i]
  test_region_index = c(1:6)[i]
  RMSE <- getRMSEforCitys(train_region_index = train_region_index, 
                          test_region_index = test_region_index
                          , only.120 = only.120
                          , ratio = ratio)
  return(RMSE)
  cat("nrow",row,"train:",train_region_index,"|| test:",test_region_index," || RMSE:", RMSE[row])
}


print(Sys.time()); start_time = Sys.time() # output: a (len-rep) list of len-n_method lists. If Pareto, then a len_n_weight list of RC = matrix(NA, 2, n_city)
cl <- parallel::makeCluster(detectCores(), outfile = "")
doParallel::registerDoParallel(cl)
deviances = foreach(i = 1:6, .combine = "c", .packages = c("ggplot2")) %dopar% {
  leaveOneCV(i
             , only.120 = T
             , ratio = T )
}
print(Sys.time() - start_time)
stopCluster(cl)

# Output average error ratio in C2
average_error_ratio = c(average_error_ratio=mean(abs(deviances)) )
saveRDS(average_error_ratio, "./output/C2_err.rds")
average_error_ratio = readRDS("./output/C2_err.rds");average_error_ratio
# SECTION 4.2 -------------------------------------------------------------


# INitialize parameters
t0_train = 12; t0_test = 12
alpha = .01;  t_end = 120 

train_region_index = test_region_index = region_index

getdeviationforCitys <- function(train_region_index, test_region_index,
                                 t0_train=12, alpha = .1
                                 , t_end = 120 
                                 , return.res = "ratio"
                                 ){
  theta = getSIRMAP(data=data, region_index = train_region_index, t0 = t0_train)
  cat("R0:", theta[13:15],"\n")
  deviations <- c()
  t0_test = t0_train
  for (l in test_region_index){
    p = plotNewValidation(data=data, l=l, t0=t0_test, theta=theta, t_end= t_end
                          , alpha = alpha, lag=9, date_break="30 day")
    if (return.res == "last") {
      deviations = c(deviations, p$deviations[length(p$deviations)] )
    }else if (return.res == "ratio") {
      deviations = c(deviations, p$deviations[length(p$deviations)] / p$dataobs$count[length(p$dataobs$count)] )
    }else{
      dev.l = c(p$dataest$mean[length(p$dataest$mean)], p$dataobs$count[length(p$dataobs$count)] ) 
      deviations = c(deviations, dev.l)
    }
  }
  return(deviations)
}

a1 = getdeviationforCitys(train_region_index=region_index, test_region_index=region_index)
a1


# return(mean(abs(deviations)))
