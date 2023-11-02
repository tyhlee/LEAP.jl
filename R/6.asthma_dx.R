library(tidyverse)

# required: pop, life_table, asthma inc, asthma prev
# goal: estimate the "un"diagnosis prob, 
# i.e., the prob that an asthma patient "outgrows" asthma and 
# become undiagnosed with asthma

#iterative method
chosen_province <- "CA"
starting_year <- 2000
end_year <- 2065

asthma_tuner <- function(chosen_province="CA",
                         starting_year= 1999,
                         end_year = 2065){
  life_table <- read_csv("../src/processed_data/master_life_table.csv") %>% 
    filter(province==chosen_province) %>% 
    filter(year >= starting_year) %>% 
    select(year,age,prob_death,sex) %>% 
    pivot_wider(names_from=sex,values_from=prob_death) %>% 
    group_split(year)
  
  asthma_prev <- read_csv("master_asthma_prev_interpolated.csv") %>% 
    filter(province==chosen_province) %>% 
    filter(year>=starting_year) %>% 
    select(-province) %>% 
    # select(year,age,pred,sex) %>% 
    # mutate(sex=ifelse(sex==1,"M","F")) %>% 
    # pivot_wider(names_from=sex,values_from=pred) %>% 
    group_split(year) %>% 
    lapply(.,function(df){
      tmp_df <- expand.grid(year=df$year %>% unique(),
                            age=96:110) %>% 
        as.data.frame() %>% 
        mutate(M=df$M[nrow(df)],
               `F`=df$`F`[nrow(df)])
      
      rbind(df,tmp_df)
    })
  
  # CA asthma incidence assumed to be equal to BC
  asthma_inc <- read_csv("master_asthma_inc_interpolated.csv") %>% 
    filter(province==chosen_province) %>% 
    filter(year>=starting_year) %>% 
    select(-province) %>% 
    # select(year,age,pred,sex) %>% 
    # mutate(sex=ifelse(sex==1,"M","F")) %>% 
    # pivot_wider(names_from=sex,values_from=pred) %>% 
    group_split(year) %>% 
    lapply(.,function(df){
      tmp_df <- expand.grid(year=df$year %>% unique(),
                            age=96:110) %>% 
        as.data.frame() %>% 
        mutate(M=df$M[nrow(df)],
               `F`=df$`F`[nrow(df)])
      
      rbind(df,tmp_df)
    })
  
  max_year_asthma_data <- max(asthma_prev[[length(asthma_prev)]])
  tmp_years <- (max_year_asthma_data+1):end_year
  tmp_asthma_inc <- asthma_inc[[length(asthma_inc)]]
  tmp_asthma_prev <- asthma_prev[[length(asthma_prev)]]
  for(l in tmp_years){
    asthma_inc <- c(asthma_inc,list(tmp_asthma_inc %>% 
                                      mutate(year=l)))
    asthma_prev <- c(asthma_prev,list(tmp_asthma_prev %>% 
                                        mutate(year=l)))
  }
  
  
  # pop_est <- read_csv("../src/processed_data/initial_pop_distribution.csv") %>% 
  #   filter(province==chosen_province) %>% 
  #   filter(year >= starting_year) %>% 
  #   mutate(sex=ifelse(sex,"Male","Female"))
  
  results_assessment <- results_dx <- results_mis <- results_calibrated_inc <- c()
  years <- c((starting_year+1):end_year)
  ages <- c(4:110)
  
  for( i in 1:length(years)){
    tmp_year <- years[i]
    tmp_prev_past <- asthma_prev[[i]] %>% 
      filter(age!=110) %>% 
      select(-year,-age)
    tmp_prev <- asthma_prev[[i+1]] %>% 
      filter(age!=3) %>% 
      select(-year,-age)
    tmp_inc <- asthma_inc[[i+1]] %>% 
      filter(age!=3) %>% 
      select(-year,-age)
    
    # final eqn
    # prev = past asthma patients remain Dxed with asthma  + non-asthma patients get asthma and Dxed with asthma
    # + non-asthma patients do not get asthma but get Dxed with asthma
    # prev = prev_past*p(reassessment) + 
    # inc*(1-prev_past)*p(correct diagnosis) + 
    # (1-inc)*(1-prev_past)*p(mis diagnosis)
    # constraints:
    # 0 <= p's <= 1
    # Consequently, we need to calibrate incidence
    
    # Here is one way to estimate the parameters; this does not reflect the true misDx rates
    
    # step 1: estimate tmp_u by assuming p(correct diagnosis) = 1 
    # tmp_prev = tmp_prev_past*tmp_u + tmp_inc*(1-tmp_prev_past)  # 0<=tmp_u<=1
    tmp_u <- (tmp_prev - tmp_inc * (1-tmp_prev_past))/tmp_prev_past
    
    # key step: bound it!
    tmp_u[tmp_u>1] <- 1
    tmp_u[tmp_u<0] <- 0
    
    # step 2: estimate tmp_d
    # tmp_prev = tmp_prev_past*tmp_u + tmp_inc*(1-tmp_prev_past)*tmp_d
    tmp_d <- (tmp_prev-tmp_prev_past*tmp_u)/(tmp_inc*(1-tmp_prev_past))
    
    # again bound it!
    tmp_d[tmp_d>1] <- 1
    tmp_d[tmp_d<0] <- 0
    
    # step 3: estimate tmp_mis
    
    # tmp_prev = tmp_prev_past * tmp_u + tmp_inc*(1-tmp_prev_past)*tmp_d +  (1-tmp_inc)*(1-tmp_prev_past)*tmp_m
    tmp_mis <- ((tmp_prev - tmp_prev_past*tmp_u - tmp_inc*(1-tmp_prev_past)*tmp_d))/((1-tmp_inc)*(1-tmp_prev_past))
    
    # should not need to bound it . . .
    # sanity check
    if(any(tmp_mis>1 | tmp_mis < 0)){
      print(paste0(chosen_province," ",i))
    }
    
    # sanity check
    try(max(abs(tmp_prev-(tmp_prev_past*tmp_u + tmp_inc*(1-tmp_prev_past)*tmp_d  + (1-tmp_inc)*(1-tmp_prev_past)*tmp_mis )))>
          1e-10, stop(paste0(i,"wrong wrong wrong")))
    
    # tmp_prev : 0.0798, 0.0371
    # tmp_prev_past : 0.0577, 0.0248
    # tmp_inc :
    # 0.105, 0.0621
    # 0.0313, 0.0819
    # 0.120, 0.0724
    # 0.895 non asthma => 0.0280135 asthma
    # 0.105 asthma =>
    # 0.120 = 0.0280135 + 0.105*0.876
    #
    # j <- 5
    #
    # tmp_prev_past[j,1]
    # tmp_inc[j,1]
    #
    # (1-0.141)*0.016 + 0.141
    # tmp_prev[j,1]
    
    tmp_df <- data.frame(year=tmp_year,
                         age=ages)
    tmp_df_2 <- data.frame(year=tmp_year,
                           age=c(3,ages))
    results_assessment[[i]] <- cbind(tmp_df,tmp_u)
    results_dx[[i]] <- cbind(tmp_df,tmp_d)
    results_mis[[i]] <- cbind(tmp_df,tmp_mis)
    
    # inc*(1-prev_past)*p(correct diagnosis) + 
    # (1-inc)*(1-prev_past)*p(mis diagnosis)
    
    results_calibrated_inc[[i]] <- cbind(tmp_df_2,rbind(asthma_inc[[i+1]] %>% 
                                                          filter(age==3) %>% 
                                                          select(-year,-age),tmp_inc*(1-tmp_prev_past)*tmp_d + 
                                                          (1-tmp_inc)*(1-tmp_prev_past)*tmp_mis))
  }
  return(list(results_assessment,results_dx,results_mis,results_calibrated_inc))
}


CA_tuner <- asthma_tuner()
BC_tuner <- asthma_tuner(chosen_province = "BC",end_year=2043)

master_assessment <- rbind(do.call(rbind,CA_tuner[[1]]) %>% mutate(province='CA'),
                           do.call(rbind,BC_tuner[[1]]) %>% mutate(province='BC'))

master_dx <- rbind(do.call(rbind,CA_tuner[[2]]) %>% mutate(province='CA'),
                   do.call(rbind,BC_tuner[[2]]) %>% mutate(province='BC'))

master_mis <- rbind(do.call(rbind,CA_tuner[[3]]) %>% mutate(province='CA'),
                    do.call(rbind,BC_tuner[[3]]) %>% mutate(province='BC'))

master_calibrated_inc <- rbind(do.call(rbind,CA_tuner[[4]]) %>% mutate(province='CA'),
                               do.call(rbind,BC_tuner[[4]]) %>% mutate(province='BC')) %>% 
  as.data.frame()

master_assessment <- master_assessment %>% select(year,age,`F`,M,province)
master_dx <- master_dx %>% select(year,age,`F`,M,province)
master_mis <- master_mis %>% select(year,age,`F`,M,province)

write_csv(master_assessment,"master_asthma_assessment.csv")
write_csv(master_dx,"master_asthma_dx.csv")
write_csv(master_mis,"master_asthma_mis_dx.csv")
# write_csv(master_calibrated_inc,"master_asthma_inc_calibrated.csv")
