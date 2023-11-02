library(tidyverse)
library(epitools)
library(mgcv)

logit <- function(p){
  log(p/(1-p))
}

inverse_logit <-function(x){
  exp(x)/(1+exp(x))
}

prev_calibrator <- function(target_prev,
                            target_OR,
                            p_risk,
                            beta0 = NULL,
                            multiple_risk_factors=F,
                            chosen_trace=F){
  if(is.null(beta0)){
    beta0 <- logit(target_prev)
  }
  
  obj_function <- function(x){
    # binary
    if(!multiple_risk_factors){
      if(length(target_OR)==1){
        p0 <- inverse_logit(beta0-p_risk*x)
        p0star <- inverse_logit(logit(p0) + log(target_OR))
        a <- p_risk * p0star
        c <- (1-p_risk) * p0
        return(abs(a + c - target_prev))
      } 
      else{
        logit_p0 <- beta0 - sum(p_risk[-1]*x)
        px <- inverse_logit(logit_p0 + log(target_OR))
        
        return(abs(sum(px * p_risk) - target_prev))
      }}
    else{
      # number of risk factors
      num_p <- length(p_risk)
      # number of levels of risk factors
      p_length <- lapply(p_risk,length) %>% unlist()
      p_length_optim <- p_length - 1 
      
      # break up x
      indices <- xs <- c()
      index_start <- 1
      for(i in 1:num_p){
        index_end <- index_start + p_length_optim[[i]] - 1 
        xs[[i]] <- x[index_start:index_end]
        index_start <- index_end+1
      } 
      
      mapply(function(tmp_p_risk,tmp_x){
        sum(tmp_p_risk[-1]*tmp_x)
      },p_risk,xs,SIMPLIFY = F) %>% unlist() -> penalty
      
      logit_p0 <- beta0 - sum(penalty)
      
      px_unlisted <- inverse_logit(logit_p0 + log(unlist(target_OR)))
      p_risk_unlisted <- unlist(p_risk)
      return(abs(sum(px_unlisted * p_risk_unlisted) - target_prev))
    }
    
    
  }
  if(!multiple_risk_factors){
    num_params <- ifelse(length(target_OR)==1,1,length(target_OR)-1)
  } else{
    num_params <- sum(unlist(target_OR)!=1)
  }
  return(optim(rep(0,num_params),obj_function,
               control=list(abstol=1e-15,maxit=10000,trace=chosen_trace),
               method="BFGS",hessian = T)$par) 
}

# ra, misDX, and Dx should be chosen ahead of time such that
# prev = past_prev*ra + inc*(1-past_prev)*Dx + (1-inc)*(1-past_prev)*misDx)
inc_calibrator <- function(target_inc,
                           past_target_prev,
                           past_target_OR,
                           target_OR,
                           p_risk,
                           ra=1,
                           misDx=0,
                           Dx=1,
                           chosen_trace=F,
                           risk_set,
                           initial_param_values=NULL,
                           method="nlm"){
  
  beta0 <- logit(target_inc)
  
  prop_asthma <- past_target_prev # (ref_b0+ref_d0)
  prop_no_asthma <- (1-past_target_prev) # (ref_a0+ref_c0)
  
  # reconstruct contingency table for each OR
  ref_p_risk <- p_risk[1]
  sol <- prev_calibrator(past_target_prev,past_target_OR, p_risk)
  p0 <- inverse_logit(logit(past_target_prev) - sum(p_risk[-1]*sol))
  calibrated_p <- inverse_logit(logit(p0) + log(past_target_OR))
  # distribution of the risk factors for the population without asthma
  no_asthma_p_risk_dist <- (1-calibrated_p) * p_risk
  # normalize
  no_asthma_p_risk_dist <- no_asthma_p_risk_dist/sum(no_asthma_p_risk_dist)
  
  # for each OR, we need to obtain the contingency table
  
  prev_table <- c()
  for(i in 1:(length(past_target_OR)-1)){
    # print(i)
    tmp_p_risk <- p_risk[c(1,i+1)]
    tmp_p_risk <- tmp_p_risk/sum(tmp_p_risk)
    tmp_p <- calibrated_p[c(1,i+1)]
    # return: a b c d
    # a: no exp, no asthma
    # b: no exp, yes asthma
    # c: yes exp, no asthma
    # d: yes exp, yes asthma
    # Solve the following:
    # tmp_target_prev  = (b+d)/(a+b+c+d) 
    # tmp_p[1] = b/(a+b) 
    # tmp_p[2] = d/(c+d) 
    # tmp_target_OR  = (a*d)/(b*c) 
    
    # someone else has done it;
    # use the metafor pkg
    # Bonett, D. G. (2007).
    # Transforming odds ratios into correlations for meta-analytic research. 
    # American Psychologist, 62(3), 254–255. ⁠https://doi.org/10.1037/0003-066x.62.3.254⁠
    
    nn <- 1e10
    prev_table[[i]] <-  rev(metafor::conv.2x2(ori=past_target_OR[i+1],
                                              ni = nn,
                                              # prev of exposure
                                              n1i = ((1-tmp_p[2])*tmp_p_risk[2] +tmp_p_risk[2] * tmp_p[2])*nn,
                                              # prev of asthma
                                              n2i=  sum(tmp_p_risk * tmp_p)*nn)/nn)
  }
  
  future_prev_table <- c()
  target_prev <- (past_target_prev*ra + target_inc*(1-past_target_prev)*Dx +
                    (1-target_inc)*(1-past_target_prev)*misDx)
  tmp_sol <- prev_calibrator(target_prev,target_OR, p_risk)
  tmp_p0 <- inverse_logit(logit(target_prev) - sum(p_risk[-1]*tmp_sol))
  tmp_calibrated_p <- inverse_logit(logit(tmp_p0) + log(target_OR))
  
  for(i in 1:(length(past_target_OR)-1)){
    # print(i)
    tmp_p_risk <- p_risk[c(1,i+1)]
    tmp_p_risk <- tmp_p_risk/sum(tmp_p_risk)
    tmp_p <- tmp_calibrated_p[c(1,i+1)]
    
    nn <- 1e10
    future_prev_table[[i]] <-  rev(metafor::conv.2x2(ori=target_OR[i+1],
                                                     ni = nn,
                                                     # prev of exposure
                                                     n1i = ((1-tmp_p[2])*tmp_p_risk[2] +tmp_p_risk[2] * tmp_p[2])*nn,
                                                     # prev of asthma
                                                     n2i=  sum(tmp_p_risk * tmp_p)*nn)/nn)
  }
  
  obj_function <- function(y){
    # break up params
    ncol_risk_set <- ncol(risk_set)
    tmp_risk_set <- risk_set
    if(ncol_risk_set!=1){
      # break up xs
      levels <- apply(risk_set,2,function(x){ length(unique(x))-1})
      tmp_x <- c()
      start_index <- 1
      for(j in 1:ncol_risk_set){
        end_index <- start_index + levels[j]-1
        tmp_x[[j]] <- risk_set %>% 
          select(colnames(risk_set)[j]) %>% 
          distinct() %>% 
          mutate(log_OR = c(0,y[start_index:end_index]))
        start_index <-  start_index + levels[j]
        tmp_risk_set <- tmp_risk_set %>% 
          left_join(tmp_x[[j]] ,by=colnames(risk_set)[j])
      }
      x <- rowSums(tmp_risk_set[,-c(1:ncol_risk_set)])[-1]
      
    } else{
      x <- y
    }
    
    # x = log(OR) for incidence eqn
    q <- length(no_asthma_p_risk_dist)-1
    target_OR_no_ref <- target_OR[-1]
    # # calibrate the current inc to the target inc
    tmp_sol <- prev_calibrator(target_inc,target_OR = exp(c(0,x)),p_risk = no_asthma_p_risk_dist)
    logit_p0 <- beta0 - sum(tmp_sol*no_asthma_p_risk_dist[-1])
    # logit_p0 <- beta0
    calibrated_inc <- inverse_logit(logit_p0 + c(0,x))
    # the following eqn should be satisfied automatically by construction
    # target_prev = (past_target_prev*ra + target_inc*(1-past_target_prev)*Dx +
    #       (1-target_inc)*(1-past_target_prev)*misDx)
    
    ref_cal_inc <- calibrated_inc[1]
    calibrated_inc_no_ref <- calibrated_inc[-1]
    
    result <- 0
    
    for(i in 1:(length(target_OR)-1)){
      cal_inc <- calibrated_inc_no_ref[i]
      target_x <- x[i]
      
      ref_a0 <- prev_table[[i]][1]
      ref_b0 <- prev_table[[i]][2]
      ref_c0 <- prev_table[[i]][3]
      ref_d0 <- prev_table[[i]][4]
      
      # contingency table of the population with asthma from a previous year
      # if ra=1, no reversibility
      a0 <- ref_b0*(1-ra)
      c0 <- ref_d0*(1-ra)
      b0 <- ref_b0*ra
      d0 <- ref_d0*ra
      
      # contingency table of the exposure level 
      # no exposure & no asthma: did not get asthma and did not get misdx + got asthma but misDx
      a1 <-  (1-ref_cal_inc)*ref_a0*(1-misDx) + ref_cal_inc * ref_a0 *(1-Dx)
      # no exposure & yes asthma: get asthma and correctly Dx + did not get asthma but misDx
      b1 <- (1-ref_cal_inc)*ref_a0*misDx + ref_cal_inc * ref_a0 * Dx
      # yes exposure & no asthma: got asthma but incorrectly Dx + did not get asthma and correctly Dx
      c1 <- (1-cal_inc)*ref_c0*(1-misDx) + cal_inc * ref_c0 * (1-Dx)
      # yes exposure & yes asthma: got asthma and correctly Dx
      d1 <-   (1-cal_inc)*ref_c0*misDx + cal_inc* ref_c0 * Dx
      # two targets
      #  objective: asthma prev OR
      a <- a0 + a1
      b <- b0 + b1
      c <- c0 + c1
      d <- d0 + d1
      tmp_OR <- a*d/(b*c)
      result <- result +
        abs(log(target_OR_no_ref[i]) - (log(d) + log(a) - log(b) - log(c)))
        # sum(abs(c(a,b,c,d)-future_prev_table[[i]]))
      
      # abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1)) #figure out why this works
      # abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1)) #figure out why this works
      
      # result <- result + abs(log(target_OR_no_ref[i]) + log((b0+b1)/b1) + log(d1/(d1+d0)) + log(a1/(a0+a1)) + log((c0+c1)/c1) - target_x)
    }
    
    # print(tmp_sol)
    return(result %>% unlist())
  }
  
  n_par <-   sum(risk_set %>% 
                   apply(.,2,function(x){length(unique(x))-1}))
  
  if(is.null(initial_param_values)){
    param_values <- rep(0,n_par)
  } else{
    param_values <- initial_param_values
  }

  if(method=="nlm"){
  tmp_sol1 <- nlm(obj_function,
                  param_values,
                  print.level=ifelse(chosen_trace,1,0),
                  iterlim = 1e5,
                  steptol=1e-19)
  tmp_sol1_par <- tmp_sol1$estimate
  } else if (method=="BFGS"){
    tmp_sol1 <- optim(param_values,
                      obj_function,
                      control=list(maxit=1e4,trace=chosen_trace,
                                   ndeps=rep(1e-10,n_par)),
                      method="BFGS")
    tmp_sol1_par <- tmp_sol1$par
  } else if (method=="L-BFGS-B"){
    tmp_sol1 <- optim(param_values,
                      obj_function,
                      control=list(trace=chosen_trace,ndeps=rep(1e-8,n_par)),
                      lower=rep(1e-10,n_par),
                      method="L-BFGS-B")
    tmp_sol1_par <- tmp_sol1$par
  }
  
  ncol_risk_set <- ncol(risk_set)
  tmp_risk_set <- risk_set
  y <- tmp_sol1_par
  
  if(ncol_risk_set!=1){
    # break up xs
    levels <- apply(risk_set,2,function(x){ length(unique(x))-1})
    tmp_x <- c()
    start_index <- 1
    for(j in 1:ncol_risk_set){
      end_index <- start_index + levels[j]-1
      tmp_x[[j]] <- risk_set %>% 
        select(colnames(risk_set)[j]) %>% 
        distinct() %>% 
        mutate(log_OR = c(0,y[start_index:end_index]))
      start_index <-  start_index + levels[j]
      tmp_risk_set <- tmp_risk_set %>% 
        left_join(tmp_x[[j]] ,by=colnames(risk_set)[j])
    }
    x <- rowSums(tmp_risk_set[,-c(1:ncol_risk_set)])[-1]
    
  } else{
    x <- y
  }
  
  
  # adjust for incidence rate
  tmp_sol2 <- prev_calibrator(target_inc,c(1,exp(x)),no_asthma_p_risk_dist,
                              beta0 = logit(target_inc) ,chosen_trace=F)
  
  return(list(tmp_sol1,tmp_sol2,
              corrector=sum(no_asthma_p_risk_dist[-1]*tmp_sol2),
              inc_OR = c(1,exp(x)),
              calibrated_inc=inverse_logit(logit(target_inc) - sum(no_asthma_p_risk_dist[-1]*tmp_sol2) + c(0,x))))
}

# return 1) correction term for prev
#        2) correction term for inc
#        3) OR for inc

asthma_prev_inc_calibrator <- function(year,age,sex,
                                       risk,
                                       OR){
     return(NA)
}
