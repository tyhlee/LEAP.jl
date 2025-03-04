library(tidyverse)
library(here)
library(mgcv)
library(roptim)
source("calibration_helper_function.R")
options(dplyr.summarise.inform = FALSE)

# inputs
chosen_province <- "CA"

max_cal_year <- 2065 # 2065 for CA; #2043 for BC
min_cal_year <- 2000
stabilization_year <- 2025

# asthma prev and inc -----------------------------------------------------

asthma_inc_model <- read_rds("asthma_incidence_model.rds")
asthma_prev_model <- read_rds("asthma_prevalence_model.rds")
asthma_max_age <- 62

asthma_predictor <- function(age,sex,year,type){
  age <- pmin(age,asthma_max_age)
  
  year <- pmin(year,stabilization_year)
  
  if(type == "prev" ){
    return( exp(predict(asthma_prev_model,newdata=data.frame(age,sex,year))) %>% 
              unlist())
  } else{
    return( exp(predict(asthma_inc_model,newdata=data.frame(age,sex,year))) %>% 
              unlist())
  }
  
}

xs <- expand.grid(age=3:110,sex=c(0,1),year=min_cal_year:max_cal_year) %>% 
  as.data.frame()

df_asthma <- xs %>% 
  mutate(inc = asthma_predictor(age,sex,year,"inc")) %>% 
  mutate(prev = asthma_predictor(age,sex,year,"prev")) %>% 
  mutate(inc = ifelse(age==3,prev,inc))

# inc <- read_csv("master_asthma_inc_interpolated.csv") %>% 
#   filter(province==chosen_province)
inc <- df_asthma %>% 
  select(year,age,sex,inc) %>% 
  pivot_wider(names_from=sex,values_from=inc) %>% 
  as.data.frame()
colnames(inc)[c(3,4)] <- c("F","M")
inc$province <- chosen_province

# prev <- read_csv("master_asthma_prev_interpolated.csv")%>% 
#   filter(province==chosen_province)
prev <- df_asthma %>% 
  select(year,age,sex,prev) %>% 
  pivot_wider(names_from=sex,values_from=prev) %>% 
  as.data.frame()
colnames(prev)[c(3,4)] <- c("F","M")
prev$province <- chosen_province

# misDx <- read_csv("master_asthma_mis_dx.csv")%>% 
#   filter(province==chosen_province)
# misDx$`F` <- 0
# misDx$M <- 0
# Dx <- read_csv("master_asthma_dx.csv") %>% 
#   filter(province==chosen_province)
# Dx$`F` <- 1
# Dx$M <- 1

RA <- read_csv("../src/processed_data/master_asthma_reassessment.csv")%>% 
  filter(province==chosen_province)


# risk factors ------------------------------------------------------------

p_fam_distribution <- data.frame(fam_history=c(0,1),
                                 prob_fam=c(1-0.2927242,0.2927242))

# Abx exposure: 0 1 2 3 4 5+
# differs by year
Abx_count_model <- read_rds("BC_count_model.rds")
p_antibiotic_exposure <- function(chosen_year,chosen_sex){
  # 2025 for females
  # 2028 for males
  # sex: 0 female; 1 male
  # to cap it
  if( chosen_sex == 1){
    chosen_year <- min(2028-1,chosen_year)
  }  else{
    chosen_year <- min(2025-1,chosen_year)
  }
  nb_mu <-   exp(predict(Abx_count_model,
                         newdata= data.frame(sex=chosen_sex,
                                             year = chosen_year,
                                             N= 1,
                                             after2005 = as.numeric(chosen_year>2005)) %>% 
                           mutate(after2005year = after2005*year),
                         type='link'))
  nb_size <- exp(Abx_count_model$family$getTheta())
  tmp_p <- dnbinom(c(0:5),mu=nb_mu,size=nb_size)
  tmp_p[6] <- 1- sum(tmp_p[1:5])
  data.frame(abx_exposure= c(0:5),prob_abx=tmp_p)
}

# prev eqn OR
OR_fam_history <- list(c(1,1.13),
                       c(1, exp((log(1.13)+log(2.4))/2)),
                       c(1,2.4))
OR_fam_history<- data.frame(age=c(3,4,5),do.call(rbind,OR_fam_history))
colnames(OR_fam_history)[-1] <- c(0,1)
OR_fam_history <- pivot_longer(OR_fam_history,cols=-1,names_to="fam_history",values_to="OR_fam") %>% 
  mutate(fam_history = as.numeric(fam_history))




OR_abx_calculator <- function(age,dose,params = c(1.711+0.115,-0.225,0.053)){
  ifelse(dose==0,1,exp(sum(params*c(1,pmin(age,7),pmin(dose,3)))))
}

OR_fam_calculator <- function(age,fam_hist,params=c(log(1.13),(log(2.4)+log(1.13))/2)){
  ifelse(age<3 | fam_hist==0 | age > 7,1,exp(params[1] + params[2]*(pmin(age,5)-3)))
}

OR_risk_factor_calculator <- function(fam_hist,age,dose,params=list(c(log(1.13),(log(1.13)+log(2.4))/2-log(1.13)),c(1.711+0.115,-0.225,0.053))){
  ifelse(age<3,1,exp(log(OR_fam_calculator(age,fam_hist,params[[1]])) + log(OR_abx_calculator(age,dose,params[[2]]))))
}

# fam_history + \beta_age * (age-3) + dose()
# free parameters are : age

OR_Abx <- read_csv("dose_response_log_aOR.csv") 
colnames(OR_Abx) <- c("age",paste0("OR",c(1:5)))
OR_Abx$OR0 <- 0
OR_Abx <- OR_Abx %>% 
  select(age,OR0,OR1:OR5) %>% 
  mutate(across(contains("OR"),exp)) %>% 
  filter(age>=3) %>% 
  rbind(data.frame(age=8,OR0=1,OR1=1,OR2=1,OR3=1,OR4=1,OR5=1))

OR_Abx <- pivot_longer(OR_Abx,cols=-1,names_to="abx_exposure",values_to="OR_abx") %>% 
  mutate(abx_exposure=as.numeric(str_remove(abx_exposure,"OR")))

risk_set <- expand.grid(fam_history = c(0,1),
                        abx_exposure = c(0,1,2,3,4,5))

risk_factor_geneator <- function(chosen_year,chosen_sex,chosen_age){
  tmp_p_fam <- p_fam_distribution
  birth_year <- chosen_year - chosen_age
  tmp_abx_exposure <- p_antibiotic_exposure(max(birth_year,2000),chosen_sex)
  tmp_OR_fam <- OR_fam_history %>% 
    filter(age == min(chosen_age,5)) %>% 
    select(-age)
  tmp_OR_abx <- OR_Abx %>% 
    filter(age == min(chosen_age,8)) %>% 
    select(-age)
  # tmp_OR_abx$OR_abx[4] <- exp(sum(log(tmp_OR_abx$OR_abx)[4:6]*(tmp_abx_exposure$prob_abx[4:6]/sum(tmp_abx_exposure$prob_abx[4:6]))))
  tmp_abx_exposure$prob_abx[4] <- sum(tmp_abx_exposure$prob_abx[4:6])
  
  tmp_OR_abx <- tmp_OR_abx %>% 
    filter(abx_exposure<=3)
  tmp_abx_exposure <- tmp_abx_exposure %>% 
    filter(abx_exposure<=3)
  
  risk_set %>%
    mutate(year = chosen_year,
           sex = chosen_sex,
           age = chosen_age) %>%
    filter(abx_exposure<=3) %>% 
    left_join(tmp_p_fam,by=c("fam_history")) %>%
    left_join(tmp_abx_exposure,by=c("abx_exposure")) %>%
    left_join(tmp_OR_fam, by = c("fam_history")) %>%
    left_join(tmp_OR_abx, by =c("abx_exposure")) %>% 
    mutate(prob = prob_fam * prob_abx,
           OR = OR_abx * OR_fam) %>%
    select(fam_history,abx_exposure,year,sex,age,prob,OR)
}

# tmp <- risk_factor_geneator(2002,1,4)
# kk <- 8
# OR_risk_factor_calculator(tmp$fam_history[kk],tmp$age[kk],tmp$abx_exposure[kk])

# algorithm for the birth cohort
tmp_inc <- inc %>% 
  select(-province) %>%
  pivot_longer(3:4,values_to="inc",names_to='sex') %>% 
  mutate(sex = as.numeric(sex=="M"))

tmp_prev <- prev %>% 
  select(-province)%>% 
  pivot_longer(3:4,values_to="prev",names_to='sex')%>% 
  mutate(sex = as.numeric(sex=="M"))

tmp_RA <- RA %>% 
  select(-province)%>% 
  pivot_longer(3:4,values_to="ra",names_to='sex')%>% 
  mutate(sex = as.numeric(sex=="M"))

# tmp_misDx <- misDx %>% 
#   select(-province)%>% 
#   pivot_longer(3:4,values_to="misDx",names_to='sex')%>% 
#   mutate(sex = as.numeric(sex=="M"))
# 
# tmp_Dx <- Dx %>% 
#   select(-province)%>% 
#   pivot_longer(3:4,values_to="Dx",names_to='sex')%>% 
#   mutate(sex = as.numeric(sex=="M"))

# for each year, sex, age
# given the effects of risk factors in the incidence equation,
# spit out the loss function 
# .5989652 -0.3574636
calibrator <- function(inc_beta_parameters=c(0.3766256,-0.225),
                       chosen_year,
                       chosen_sex,
                       chosen_age){
  
  if(!is.list(inc_beta_parameters)){
    inc_beta_parameters <- list(c(log(1.13),inc_beta_parameters[1]),c(1.826,inc_beta_parameters[2],0.053))
  }
  
  if(chosen_age<=7){
    tmp_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age)
    
    target_prev <- tmp_prev %>% 
      filter(age == chosen_age & 
               year == chosen_year &
               sex == chosen_sex) %>% 
      select(prev) %>% 
      unlist()
    
    target_OR <- tmp_risk_set$OR
    target_risk_p <- tmp_risk_set$prob
    prev_sol <- prev_calibrator(target_prev,target_OR,target_risk_p)
    p0 <- inverse_logit(logit(target_prev) - sum(target_risk_p[-1]*prev_sol))
    calibrated_prev <- inverse_logit(logit(p0) + log(target_OR))
    
    tmp_risk_set$calibrated_prev <- calibrated_prev
    tmp_risk_set$prev <- target_prev
    
    if(chosen_year== 2000){
      return(tmp_risk_set)
    }
    
    if(chosen_age==3){
      tmp_risk_set$calibrated_inc <- calibrated_prev
    } 
    
    else{ # aged 4 or more
      
      target_inc <- tmp_inc %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(inc) %>% 
        unlist()
      
      tmp_risk_set$inc <- target_inc
      
      past_target_prev <- tmp_prev %>% 
        filter(age == chosen_age-1 & 
                 year == max(min_cal_year,chosen_year-1) &
                 sex == chosen_sex) %>% 
        select(prev) %>% 
        unlist()
      
      past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1)
      past_target_OR <- past_risk_set$OR
      past_target_risk_p <- past_risk_set$prob
      
      target_RA <- tmp_RA %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(ra) %>% 
        unlist()
      
      # target_Dx <- tmp_Dx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(Dx) %>% 
      #   unlist()
      
      target_Dx <- 1
      
      # target_misDx <- tmp_misDx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(misDx) %>% 
      #   unlist()
      
      target_misDx <- 0
      
      # if(!(abs(target_prev -(
      #   past_target_prev*target_RA + 
      #   (1-past_target_prev)*target_inc*target_Dx + 
      #   (1-past_target_prev)*(1-target_inc)*target_misDx))<1e-10)){
      #   return("Something wrong with Dx misDx RA")
      # }
      
      inc_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
        select(fam_history,abx_exposure,year,sex,age,prob)
      
      inc_risk_set$prob <- inc_risk_set$prob / sum(inc_risk_set$prob)
      
      inc_risk_set$OR <- inc_risk_set %>% 
        apply(.,1,FUN=function(x){
          OR_risk_factor_calculator(fam_hist=x[1],age=x[5],dose=x[2],params=inc_beta_parameters)
        })
      
      
      
      
      # target_inc = target_inc;
      # past_target_prev = past_target_prev;
      # past_target_OR = past_target_OR;
      # target_OR = target_OR;
      # p_risk = past_target_risk_p;
      # ra = target_RA;
      # misDx = target_misDx;
      # Dx = target_Dx;
      # risk_set=inc_risk_set;
      # log_inc_OR = log(inc_risk_set$OR)
      
      inc_sol <- inc_loss_function(target_inc = target_inc,
                                   past_target_prev = past_target_prev,
                                   past_target_OR = past_target_OR,
                                   target_OR = target_OR,
                                   p_risk = past_target_risk_p,
                                   ra = target_RA,
                                   misDx = target_misDx,
                                   Dx = target_Dx,
                                   risk_set=inc_risk_set,
                                   log_inc_OR = log(inc_risk_set$OR)[-1])
      
    } 
  }
  # age > 7 => OR =1 for all abx
  else{
    
    tmp_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
      group_by(fam_history,year,sex,age) %>% 
      summarise(prob = sum(prob),
                OR = mean(OR)) %>% 
      ungroup() 
    
    target_prev <- tmp_prev %>% 
      filter(age == chosen_age & 
               year == chosen_year &
               sex == chosen_sex) %>% 
      select(prev) %>% 
      unlist()
    target_OR <- tmp_risk_set$OR
    target_risk_p <- tmp_risk_set$prob
    prev_sol <- prev_calibrator(target_prev,target_OR,target_risk_p)
    p0 <- inverse_logit(logit(target_prev) - sum(target_risk_p[-1]*prev_sol))
    calibrated_prev <- inverse_logit(logit(p0) + log(target_OR))
    
    tmp_risk_set$calibrated_prev <- calibrated_prev
    tmp_risk_set$prev <- target_prev
    
    if(chosen_year== 2000){
      return(tmp_risk_set)
    } else
    { 
      
      target_inc <- tmp_inc %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(inc) %>% 
        unlist()
      tmp_risk_set$inc <- target_inc
      
      past_target_prev <- tmp_prev %>% 
        filter(age == chosen_age-1 & 
                 year == max(min_cal_year,chosen_year-1) &
                 sex == chosen_sex) %>% 
        select(prev) %>% 
        unlist()
      if(chosen_age != 8){
        past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1) %>% 
          group_by(fam_history) %>% 
          summarise(prob = sum(prob),
                    OR = mean(OR))
      } else{
        
        past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1)
        
        
        ttt_target_OR <- past_risk_set$OR
        ttt_target_risk_p <- past_risk_set$prob
        ttt_prev_sol <- prev_calibrator(past_target_prev,ttt_target_OR,ttt_target_risk_p)
        ttt_p0 <- inverse_logit(logit(past_target_prev) - sum(ttt_target_risk_p[-1]*ttt_prev_sol))
        ttt_calibrated_prev <- inverse_logit(logit(ttt_p0) + log(ttt_target_OR))
        past_risk_set$calibrated_prev <- ttt_calibrated_prev
        past_risk_set %>% 
          mutate(yes_asthma = calibrated_prev * prob,
                 no_asthma = (1-calibrated_prev) * prob) -> tmp_look
        past_tmp_OR <- sum(tmp_look$no_asthma[tmp_look$fam_history==0])*sum(tmp_look$yes_asthma[tmp_look$fam_history==1])/
          (sum(tmp_look$yes_asthma[tmp_look$fam_history==0])*sum(tmp_look$no_asthma[tmp_look$fam_history==1]))
        past_risk_set %>% 
          group_by(fam_history) %>% 
          summarise(prob=sum(prob)) ->past_risk_set
        past_risk_set$OR <- c(1,past_tmp_OR)
      }
      
      past_target_OR <- past_risk_set$OR
      past_target_risk_p <- past_risk_set$prob
      
      target_RA <- tmp_RA %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(ra) %>% 
        unlist()
      
      # target_Dx <- tmp_Dx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(Dx) %>% 
      #   unlist()
      
      target_Dx <- 1
      
      # target_misDx <- tmp_misDx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(misDx) %>% 
      #   unlist()
      
      target_misDx <- 0
      
      inc_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
        filter(abx_exposure==0) %>% 
        select(fam_history,abx_exposure,year,sex,age,prob)
      
      inc_risk_set$OR <- inc_risk_set %>% 
        apply(.,1,FUN=function(x){
          OR_risk_factor_calculator(fam_hist=x[1],age=x[5],dose=x[2],params=inc_beta_parameters)
        })
      
      
      
      inc_sol <- inc_loss_function(target_inc = target_inc,
                                   past_target_prev = past_target_prev,
                                   past_target_OR = past_target_OR,
                                   target_OR = target_OR,
                                   p_risk = past_target_risk_p,
                                   ra = target_RA,
                                   misDx = target_misDx,
                                   Dx = target_Dx,
                                   risk_set=inc_risk_set,
                                   log_inc_OR = log(inc_risk_set$OR)[-1])
    } 
  }
  
  
  return(inc_sol)
  
}

# inc_beta_parameters <-list(c(log(1.13),(log(2.4)-log(1.13))/2),c(1.711+0.115,-0.225,0.053))
# initial values
inc_beta_parameters <-c( (log(2.4)-log(1.13)) / 2 , -0.225)

# chosen_year <- 2002
# chosen_age <- 4
# chosen_sex <- 1
# calibrator(inc_beta_parameters,2004,1,10)


baseline_year=2001
stabilization_year=2025
max_age=63

inc_beta_solver <- function(baseline_year=2001,stabilization_year=2025,max_age=63){
  cal_years <- baseline_year:(stabilization_year+1)
  ages <- 4:max_age
  sexes <- 0:1
  covar <- expand.grid(year=cal_years,sex=sexes,age=ages) %>% 
    as.data.frame()
  
  
  
  obj <- function(inc_beta_parameters){
    apply(covar,1,FUN = function(x){
      calibrator(inc_beta_parameters,x[1],x[2],x[3])
    }) %>% mean()
    
  }
  
  res_optim <- optim(unlist(inc_beta_parameters),fn=obj,method='BFGS')
  res_nlm <- nlm(obj,unlist(inc_beta_parameters),steptol=1e-6,gradtol=1e-6,print.level=2)
  write_rds(res_optim,"res_optim.rds")
}


# incorporate the estimates of the risk factors and correction terms ----------
res_optim <- read_rds("res_optim.rds") 
optimized_inc_beta <- res_optim$par

cal_years <- (baseline_year-1):(stabilization_year+1)
ages <- 3:max_age
sexes <- 0:1
calibration_results <- expand.grid(year=cal_years,sex=sexes,age=ages) %>%
  as.data.frame()



calculate_correction <- function(inc_beta_parameters=optimized_inc_beta,
                       chosen_year,
                       chosen_sex,
                       chosen_age){
  
  tmp_res <- data.frame(year=chosen_year,
                        sex= chosen_sex,
                        age = chosen_age,
                        obj_value=NA,
                        prev_correction = NA,
                        inc_correction = NA)
  
  if(!is.list(inc_beta_parameters)){
    inc_beta_parameters <- list(c(log(1.13),inc_beta_parameters[1]),c(1.826,inc_beta_parameters[2],0.053))
  }
  
  if(chosen_age<=7){
    tmp_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age)
    
    target_prev <- tmp_prev %>% 
      filter(age == chosen_age & 
               year == chosen_year &
               sex == chosen_sex) %>% 
      select(prev) %>% 
      unlist()
    
    target_OR <- tmp_risk_set$OR
    target_risk_p <- tmp_risk_set$prob
    prev_sol <- prev_calibrator(target_prev,target_OR,target_risk_p)
    p0 <- inverse_logit(logit(target_prev) - sum(target_risk_p[-1]*prev_sol))
    calibrated_prev <- inverse_logit(logit(p0) + log(target_OR))
    
    prevalence_correction_term <- -sum(target_risk_p[-1]*prev_sol)
    tmp_res$prev_correction <- prevalence_correction_term
    
    tmp_risk_set$calibrated_prev <- calibrated_prev
    tmp_risk_set$prev <- target_prev
    
    if(chosen_year== 2000){
      return(tmp_res)
    }
    
    if(chosen_age==3){
      tmp_res$inc_correction <- prevalence_correction_term
      return(tmp_res)
      } 
    
    else{ # aged 4 or more
      
      target_inc <- tmp_inc %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(inc) %>% 
        unlist()
      
      tmp_risk_set$inc <- target_inc
      
      past_target_prev <- tmp_prev %>% 
        filter(age == chosen_age-1 & 
                 year == max(min_cal_year,chosen_year-1) &
                 sex == chosen_sex) %>% 
        select(prev) %>% 
        unlist()
      
      past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1)
      past_target_OR <- past_risk_set$OR
      past_target_risk_p <- past_risk_set$prob
      
      target_RA <- tmp_RA %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(ra) %>% 
        unlist()
      
      # target_Dx <- tmp_Dx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(Dx) %>% 
      #   unlist()
      
      target_Dx <- 1
      
      # target_misDx <- tmp_misDx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(misDx) %>% 
      #   unlist()
      
      target_misDx <- 0
      
      
      inc_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
        select(fam_history,abx_exposure,year,sex,age,prob)
      
      inc_risk_set$prob <- inc_risk_set$prob / sum(inc_risk_set$prob)
      
      inc_risk_set$OR <- inc_risk_set %>% 
        apply(.,1,FUN=function(x){
          OR_risk_factor_calculator(fam_hist=x[1],age=x[5],dose=x[2],params=inc_beta_parameters)
        })
      
      inc_sol <- inc_correction_calculator(target_inc = target_inc,
                                   past_target_prev = past_target_prev,
                                   past_target_OR = past_target_OR,
                                   target_OR = target_OR,
                                   p_risk = past_target_risk_p,
                                   ra = target_RA,
                                   misDx = target_misDx,
                                   Dx = target_Dx,
                                   risk_set=inc_risk_set,
                                   log_inc_OR = log(inc_risk_set$OR)[-1])
      
      tmp_res$obj_value <- inc_sol[1]
      tmp_res$inc_correction <- inc_sol[2]
      
      return(tmp_res)
      
    } 
  }
  # age > 7 => OR =1 for all abx
  else{
    
    tmp_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
      group_by(fam_history,year,sex,age) %>% 
      summarise(prob = sum(prob),
                OR = mean(OR)) %>% 
      ungroup() 
    
    target_prev <- tmp_prev %>% 
      filter(age == chosen_age & 
               year == chosen_year &
               sex == chosen_sex) %>% 
      select(prev) %>% 
      unlist()
    target_OR <- tmp_risk_set$OR
    target_risk_p <- tmp_risk_set$prob
    prev_sol <- prev_calibrator(target_prev,target_OR,target_risk_p)
    p0 <- inverse_logit(logit(target_prev) - sum(target_risk_p[-1]*prev_sol))
    calibrated_prev <- inverse_logit(logit(p0) + log(target_OR))
    
    tmp_res$prev_correction <- -sum(target_risk_p[-1]*prev_sol)
    
    tmp_risk_set$calibrated_prev <- calibrated_prev
    tmp_risk_set$prev <- target_prev
    
    if(chosen_year== 2000){
      return(tmp_res)
    } else
    { 
      
      target_inc <- tmp_inc %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(inc) %>% 
        unlist()
      tmp_risk_set$inc <- target_inc
      
      past_target_prev <- tmp_prev %>% 
        filter(age == chosen_age-1 & 
                 year == max(min_cal_year,chosen_year-1) &
                 sex == chosen_sex) %>% 
        select(prev) %>% 
        unlist()
      if(chosen_age != 8){
        past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1) %>% 
          group_by(fam_history) %>% 
          summarise(prob = sum(prob),
                    OR = mean(OR))
      } else{
        
        past_risk_set <- risk_factor_geneator(max(min_cal_year,chosen_year-1),chosen_sex,chosen_age-1)
        
        
        ttt_target_OR <- past_risk_set$OR
        ttt_target_risk_p <- past_risk_set$prob
        ttt_prev_sol <- prev_calibrator(past_target_prev,ttt_target_OR,ttt_target_risk_p)
        ttt_p0 <- inverse_logit(logit(past_target_prev) - sum(ttt_target_risk_p[-1]*ttt_prev_sol))
        ttt_calibrated_prev <- inverse_logit(logit(ttt_p0) + log(ttt_target_OR))
        past_risk_set$calibrated_prev <- ttt_calibrated_prev
        past_risk_set %>% 
          mutate(yes_asthma = calibrated_prev * prob,
                 no_asthma = (1-calibrated_prev) * prob) -> tmp_look
        past_tmp_OR <- sum(tmp_look$no_asthma[tmp_look$fam_history==0])*sum(tmp_look$yes_asthma[tmp_look$fam_history==1])/
          (sum(tmp_look$yes_asthma[tmp_look$fam_history==0])*sum(tmp_look$no_asthma[tmp_look$fam_history==1]))
        past_risk_set %>% 
          group_by(fam_history) %>% 
          summarise(prob=sum(prob)) ->past_risk_set
        past_risk_set$OR <- c(1,past_tmp_OR)
      }
      
      past_target_OR <- past_risk_set$OR
      past_target_risk_p <- past_risk_set$prob
      
      target_RA <- tmp_RA %>% 
        filter(age == chosen_age & 
                 year == chosen_year &
                 sex == chosen_sex) %>% 
        select(ra) %>% 
        unlist()
      
      # target_Dx <- tmp_Dx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(Dx) %>% 
      #   unlist()
      
      target_Dx <- 1
      
      # target_misDx <- tmp_misDx %>% 
      #   filter(age == chosen_age & 
      #            year == chosen_year &
      #            sex == chosen_sex) %>% 
      #   select(misDx) %>% 
      #   unlist()
      
      target_misDx <- 0
      
      inc_risk_set <- risk_factor_geneator(chosen_year,chosen_sex,chosen_age) %>% 
        filter(abx_exposure==0) %>% 
        select(fam_history,abx_exposure,year,sex,age,prob)
      
      inc_risk_set$OR <- inc_risk_set %>% 
        apply(.,1,FUN=function(x){
          OR_risk_factor_calculator(fam_hist=x[1],age=x[5],dose=x[2],params=inc_beta_parameters)
        })
      
      

      inc_sol <- inc_correction_calculator(target_inc = target_inc,
                                           past_target_prev = past_target_prev,
                                           past_target_OR = past_target_OR,
                                           target_OR = target_OR,
                                           p_risk = past_target_risk_p,
                                           ra = target_RA,
                                           misDx = target_misDx,
                                           Dx = target_Dx,
                                           risk_set=inc_risk_set,
                                           log_inc_OR = log(inc_risk_set$OR)[-1])
      
      tmp_res$obj_value <- inc_sol[1]
      tmp_res$inc_correction <- inc_sol[2]
      
      return(tmp_res)
    } 
  }
  
}


generate_correction <- function(tmp_df){
    apply(tmp_df,1,FUN=function(x){
      calculate_correction(inc_beta_parameters=optimized_inc_beta,
                           chosen_year=x[1],
                           chosen_sex = x[2],
                           chosen_age = x[3])
    }) %>% 
    do.call(rbind,.)
}

df_correct <- generate_correction(calibration_results)

# write_rds(df_correct,"df_correct_25.rds")

df_correct_prev <- df_correct %>% 
  select(1:3,5) %>% 
  rename(correction=prev_correction) %>% 
  mutate(type='prev')

df_correct_inc <- df_correct %>% 
  select(1:3,6) %>% 
  rename(correction=inc_correction) %>% 
  mutate(type='inc')

master_correct <- rbind(df_correct_prev,
                        df_correct_inc)

write_csv(master_correct %>% 
            mutate(correction=ifelse(is.na(correction),0,correction)),"../src/processed_data/master_asthma_occurrence_correction.csv")

# examine_results <- df_correct %>% 
#   filter(age !=3)
# 
# ggplot(examine_results %>%
#          filter(age <= 10) %>% 
#          mutate(year=as.factor(year)),aes(y=obj_value,x=age,col=year)) +
#   geom_point() + 
#   facet_grid(.~sex)
# # write_csv(final_result,"master_calibrated_asthma_prev_inc_M3.csv")
# 
# #### post-processing
# final_result <- read_csv("master_calibrated_asthma_prev_inc_M3.csv")
# 
# tmp <- head(final_result %>% filter(year==2004),n=8)
# tmp$calibrated_prev[seq(2,8,by=2)]/tmp$calibrated_prev[seq(1,7,by=2)]
# 
# logit <- function(p){
#   log(p/(1-p))
# }
# 
# inv_logit <- function(x){
#   exp(x)/(1+exp(x))
# }
# 
# final_result %>% 
#   group_by(fam_history,abx_exposure,sex,age) %>% 
#   summarise(mean(OR),median(OR),sd(OR)) %>% 
#   filter(`sd(OR)`!=0)  # variations are negligible 
# 
# tmp <- final_result %>% 
#   filter(age==5 & sex==0 & fam_history==0 & abx_exposure==0)
# 
