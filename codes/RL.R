#### Main Algorithms ##################################

#' Algorithm 1: SPSA
#' at time t0, based on the  current estimated parameters, for city l, solve the optimal policy by planning.
#'
#' @param theta (gamma, beta)
#' @param XR_t0_lag {X^R_t}_{t0-8}^{t0}
#' @param A_t0_1_lag {A_t}_{t0-8}^{t0-1}
#' @param GDP_daily in trillion now
#' @param w weight, R^+
#' @param t0,T_end date of decision and ending dete.(only T_end - t0 matters)
#' @param M num of replications for each parameter
#' @param range1 upper searching bound for the first grid search. range1[1] for lam2 and range1[2] for lam3
#' @param range2 searching radius for the second grid search.
#' @Output:
#'    list(lam = vlam, A = A)

policySearch <- function(theta = NA, t0 = 12, T_end = 120,
                         XR_t0_lag = NA, 
                         A_t0_1_lag = NA, 
                         lag = 9, 
                         popu = 16330000,
                         GDP_daily = 796.3978 / 365,
                         range1 = c(300, 1000), step1 = c(50, 200),
                         range2 = c(50, 200), step2 = c(5, 20),
                         range3 = c(5, 20), step3 = c(1, 4),
                         range4 = c(1, 4), step4 = c(0.25, 1),
                         step5 = c(0.01, 0.2), 
                         w = 1, J = 3, 
                         M = 50, 
                         seed = 1, 
                         is.state.oracle = F, 
                         S_t0 = NA, 
                         classifer.echo = T, 
                         echo = F, early_stop_tol = 20
                         , decision.t = 7
                         , method = "grid", model = NA, MC = F, inpre = F, only_action = F
                         , catious_c = 0.7
){
  ## Initialization ----
  
  start_time = Sys.time()
  
  if (is.state.oracle) {
    belief_state_t0 = S_t0
  }else{
    belief_state_t0 = generateBeliefState(theta, XR_t0_lag, A_t0_1_lag, popu) # S, I, R
  }
  
  if (method %in% c("grid", "both")) {
    
    gridSearch <-function(lower_para1, upper_para1, step_para1, upper_para2, step_para2, lower_para2, tol = early_stop_tol){#lower_para2
      #' @param lower_para1,upper_para1: lower and upper bound for search space of parameter 1
      #' @param upper_para2,lower_para2: lower and upper bound for search space of parameter 2
      #' @param step_para1,step_para2: step sizes
      #' @details missing(lower_para2): 
      #' @details first round: [[1, upper], [vlam1, upper]] with [step1[1],step2[1]]
      #' @details second round: [[max(vlam1[1] - range2[1], 1), vlam1[1] + range2[1]], [vlam1, last_res_for_para_2 + range2[2]]]
      #' @details third round: [[last - 1, last +1], [last - 1, last + 1]] with [0.1, 0.1]
      #' 
      
      s_seed = seed * 123
      
      vlam_best = rep(0, 0)
      value_best = rep(-1e10, 3)
      count_no_increase = 0
      for (vlam1 in seq(lower_para1, upper_para1, step_para1)){
        count = 0
        value_old  = 0
        
        if (missing(lower_para2)){
          lower_para2 = vlam1
          #lower_para2 = upper_para2
        }
        for (vlam2 in seq(max(lower_para2,vlam1), upper_para2, step_para2)) {
          vlam = c(vlam1, vlam2)
          r  = getValueforThisLambda(vlam, XR_t0_lag, A_t0_1_lag, popu, t0, T_end, theta, GDP_daily, w, M, MC = MC, 
                                     decision.t = decision.t
                                     , init_seed = s_seed, Y_smallest = value_best[1]); s_seed = s_seed + 1
          value = r[1:3]
          YRC_std = r[4:6]
          if(echo)cat("vlam = ", vlam,": values (Y,R,C) = ", round(value, 2), "std in", M ,"times:", round(YRC_std, 2), "\n")
          
          if (value[1] > value_best[1]) {
            vlam_best = vlam
            value_best = value
            count_no_increase = 0
          }else{
            count_no_increase = count_no_increase + 1
          }
          # if(count_no_increase > tol)break # [think]
          
          if(value[1] == value_old[1]){
            count = count + 1
          }
          if(count == 5)break
          value_old = value
        }
        if(echo)cat("\n")
        if(count_no_increase > tol)break
      }
      
      return(list(vlam_best = vlam_best, value_best = value_best))
    }
    # if (missing(lower_para2)){lower_para2 = vlam1}
    
    vlam1 = gridSearch(lower_para1 = 0, upper_para1 = range1[1], step_para1 = step1[1], 
                       upper_para2 = range1[2], step_para2 = step1[2])$vlam_best
    if(echo)cat("round 2 from", vlam1, "\n", "----------- FINER BEGINS! ----------- ","\n")
    
    if (only_action ) {
      A = isOutRange(XI = belief_state_t0[2], 
                     lower.search.lam1 = max(vlam1[1] - range2[1], 0), upper.search.lam1 = vlam1[1] + range2[1], 
                     lower.search.lam2 = vlam1[2] - range2[2], 
                     upper.search.lam2 = vlam1[2] + range2[2])
      if (A) {
        # if (t0 %% decision.t == t0 %% decision.t)cat(t0, "->", A, "||\t") 
        return(A = A)
      }
    }
    
    # step2[1]
    vlam2 = gridSearch(lower_para1 = max(vlam1[1] - range2[1], 0), upper_para1 = vlam1[1] + range2[1], step_para1 = step2[1],
                       upper_para2 = vlam1[2] + range2[2], lower_para2 = vlam1[2] - range2[2], step_para2 = step2[2])$vlam_best
    if(echo)cat("round 3 from", vlam2, "\n", "----------- LAST SEARCH BEGINS! ----------- ","\n")
    if (only_action ) {
      A = isOutRange(XI = belief_state_t0[2], lower.search.lam1 = max(vlam2[1] - range3[1], 0), 
                     upper.search.lam1 = vlam2[1] + range3[1], 
                     lower.search.lam2 = vlam2[2] - range3[2], 
                     upper.search.lam2 = vlam2[2] + range3[2])
      if (A) {
        # if (t0 %% decision.t == t0 %% decision.t)cat(t0, "->", A, "||\t") 
        return(A = A)
      }
    }
    # step3[1]
    vlam3 = gridSearch(lower_para1 = max(vlam2[1] - range3[1], 0), upper_para1 = vlam2[1] + range3[1], step_para1 = step3[1],
                       upper_para2 = vlam2[2] + range3[2], lower_para2 = vlam2[2] - range3[2], step_para2 = step3[2])$vlam_best
    if(echo)cat("round 4 from", vlam3, "\n", "----------- LAST SEARCH BEGINS! ----------- ","\n")
    if (only_action ) {
      A = isOutRange(XI = belief_state_t0[2], lower.search.lam1 = max(vlam3[1] - range4[1], 0), 
                     upper.search.lam1 = vlam3[1] + range4[1], 
                     lower.search.lam2 = vlam3[2] - range4[2], 
                     upper.search.lam2 = vlam3[2] + range4[2])
      if (A) {
        # if (t0 %% decision.t == t0 %% decision.t)cat(t0, "->", A, "||\t") 
        return(A = A)
      }
    }
    # step4[1]
    vlam4 = gridSearch(lower_para1 = max(vlam3[1] - range4[1], 0), upper_para1 = vlam3[1] + range4[1], step_para1 = step4[1],
                       upper_para2 = vlam3[2] + range4[2], lower_para2 = vlam3[2] - range4[2], step_para2 = step4[2])$vlam_best
    if(echo)cat("round 5 from", vlam4, "\n", "----------- LAST SEARCH BEGINS! ----------- ","\n")
    res5 = gridSearch(lower_para1 = max(vlam4[1] - step4[1], 0), upper_para1 = vlam4[1] + step4[1], step_para1 = step5[1],
                      upper_para2 = vlam4[2] + step4[2], lower_para2 = vlam4[2] - step4[2], step_para2 = step5[2])
    
    vlam_best = res5$vlam_best
    value_best = res5$value_best    
    A = policy(belief_state_t0, vlam_best)
    #if (!inpre && classifer.echo) {
    # cat("vlam_best:", vlam_best,", with value_best( Y, R, C) = ", value_best, ", action =", A, "\n")
    #cat("time cost for policy search:", round(Sys.time() - start_time, 1), "\n\n")
    #}
    
  }
  
  
  if (only_action) {
    # if (t0 %% decision.t == t0 %% decision.t)cat(t0, "->", A, "||\t")
    return(A)
  }else{
    if (!inpre && classifer.echo && method == "grid") {
      cat("vlam_best:", vlam_best,", with value_best( Y, R, C) = ", value_best, ", action =", A, "\n")
      cat("time cost for policy search:", round(Sys.time() - start_time, 1), "\n\n")
    }
    return(list(lam = vlam_best, A = A))
  }
}


#' isOutRange 
#'
#'@details if smaller than the lower search bound of para1, then must be action 1
#'@details if larger than the upper search bound of para2, then must be action 3
isOutRange <- function(XI, lower.search.lam1, upper.search.lam1, lower.search.lam2, upper.search.lam2){
  if (XI < lower.search.lam1)return(1)
  if (XI > upper.search.lam2)return(3)
  if(XI > upper.search.lam1 & XI < lower.search.lam2)return(2)
  return(F)
}
