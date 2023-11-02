library(tidyverse)
library(here)
library(mgcv)
library(roptim)

# we need: 
# 1) asthma prev
# 2) target rate (this is in per pop)
# 3) population 
# 4) 2) + 3) => annual number of severe hospitalizations needed
# 5) 1) + exacerbation module => annual number of exacerbations
# 6) so we can identify a tuner to do this

# asthma prev
chosen_province <- "BC"
min_cal_year <- 2000
baseline_year <- 2000
if(chosen_province=="CA"){
max_cal_year <- 2065

} else{
  max_cal_year <- 2043
}
chosen_projection_scenario <- "M3"
growth_type <- "M3"

prev <- read_csv("master_asthma_prev_interpolated.csv") %>% 
  filter(province==chosen_province)

inc  <- read_csv("master_asthma_prev_interpolated.csv") %>% 
  filter(province==chosen_province)

tmp_prev <- prev%>% filter(year == max(prev$year))
impute_years <- (max(prev$year)+1):(max_cal_year)

for(i in impute_years){
  prev <- rbind(prev,
                tmp_prev %>% 
                  mutate(year=i))
}

prev <- prev %>% 
  select(-province) %>% 
  pivot_longer(-c(1,2),names_to="sex",values_to="prev") %>% 
  mutate(sex=as.numeric(sex=="M"))

# target population
df_cihi <- read_rds(paste0("public_dataset/asthma_hosp/",
                           chosen_province,"/tab1.rds"))$rate %>% 
  filter(fiscal_year >= baseline_year) %>% 
  rename(year=fiscal_year) %>% 
  pivot_longer(-1,names_to="type",values_to = "true_rate") %>% 
  filter(!is.na(true_rate)) %>% 
  mutate(type=str_remove(type,"\\+")) %>% 
  mutate(sex=case_when(str_detect(type,"M") ~ 1,
                       str_detect(type,"F") ~ 0,
                       TRUE ~ NA)) %>% 
  filter(!is.na(sex)) %>% 
  mutate(age=parse_number(type)) %>% 
  filter(!is.na(age)) %>% 
  select(-type) %>% 
  filter(age>=3) %>% 
  arrange(year,sex,age)

tmp_cihi <- df_cihi%>% filter(year == max(df_cihi$year))
impute_years <- (max(df_cihi$year)+1):(max_cal_year)

for(i in impute_years){
  df_cihi <- rbind(df_cihi,
                   tmp_cihi %>% 
                  mutate(year=i))
}


# population
pop_est <- read_csv("../src/processed_data/initial_pop_distribution.csv") %>% 
  filter(province==chosen_province) %>% 
  filter(year >= baseline_year) %>% 
  mutate(sex=ifelse(sex,"Male","Female"))

# grwoth type
pop_CA_BC <- read_rds(here("../src","processed_data","pop_projection_BC_CA.rds")) %>%
  select(REF_DATE,GEO,Projection_scenario,Sex,Age_group,VALUE) %>% 
  mutate(Projection_scenario = str_remove(Projection_scenario,"Projection scenario "),
         Projection_scenario = str_remove(Projection_scenario, "\\:.*")) %>% 
  filter(GEO %in% c('Canada',"British Columbia")) %>% 
  filter(REF_DATE>=(max(pop_est$year)+1)) %>% 
  rename(year=REF_DATE,
         province=GEO,
         projection_scenario = Projection_scenario,
         sex = Sex,
         ag = Age_group,
         n = VALUE) %>% 
  mutate(province = ifelse(province=="Canada","CA","BC"),
         sex  = substr(sex,1,1),
         age_remove = str_detect(ag,'to|over|All|Median|Average'),
         age_keep = ag=="100 years and over") %>% 
  filter(!age_remove | age_keep) %>%
  mutate(age = as.numeric(gsub("([0-9]+).*$","\\1",ag)),
         age = ifelse(is.na(age),0,age)) %>% 
  select(year,sex,age,province,n,projection_scenario) %>% 
  filter(sex!="B") %>% 
  mutate(n=n*1000) %>% 
  filter(province==chosen_province) %>% 
  filter(!is.na(n)) %>% 
  filter(projection_scenario == growth_type) %>% 
  select(-projection_scenario) %>% 
  mutate(sex=ifelse(sex=="M","Male","Female"))

pop <- rbind(pop_est %>% 
               select(-total_n,-prop),
             pop_CA_BC) %>% 
  filter(year <= max_cal_year) %>% 
  filter(age>=3) %>% 
  mutate(age=ifelse(age>90,90,age)) %>% 
  group_by(year,sex,age) %>% 
  summarise(n=sum(n)) %>% 
  ungroup() %>% 
  mutate(sex=as.numeric(sex=="Male"))

# control level eqn

logistic <- function(x){
  if(is.infinite(exp(x))){
    return(1)
  }
  exp(x)/(1+exp(x))
}

control_prediction <- function(year,sex,age,theta=c(-1e5,-0.3950, 2.754,1e5)){
  age_scaled = age / 100
  eta <- age_scaled*3.5430381 + sex*0.2347807 + age_scaled * sex *-0.8161495 +
    age_scaled^2 * sex *-1.1654264 + age_scaled^2 * -3.4980710
  results <- rep(0,3)
  for(j in 1:(length(theta)-1)){
    results[j] <- logistic(theta[j+1]-eta) - logistic(theta[j]-eta)
  }
  results
}

exacerbation_prediction <- function(year,sex,age,
                                    beta_control = log(c(0.1880058,
                                                         0.3760116,
                                                         0.5640174))){
  if(age<3){
    return(0)
  }
  control <- control_prediction(year,sex,age)
  exp(sum(control * beta_control))
}

p_hosp <- 0.026
# # rate multiplier given past history of hosp
# ped_RM <- 1.79
# # rate multiplier given past history of hosp
# adult_RM <- 2.88
# 
# # distribution of first asthma age
# # chance of having at least one hosp
# # avg number of total exacer
# prev_exacerbation_prediction(year,sex,age){
#   if(age==3){
#     return(0)
#   } else{
#     #
#   }
# }



df_target <- pop %>% 
  left_join(df_cihi,by=c("year","sex","age")) %>% 
  mutate(true_n = true_rate * n / 100000) %>% 
  left_join(prev,by=c("year","sex","age")) %>% 
  mutate(n_asthma = prev*n) %>% 
  rowwise() %>% 
  mutate(mean_annual_exacerbation = exacerbation_prediction(year,sex,age),
         expected_exacerbations = mean_annual_exacerbation*n_asthma,
         expected_n= p_hosp*expected_exacerbations,
         calibrator_multiplier = true_n/expected_n)

exac_cal <- df_target %>% select(year:age,calibrator_multiplier) %>% mutate(province="BC")

write_csv(exac_cal,"exacerbation_calibration_BC.csv")

final_result <- rbind(read_csv("exacerbation_calibration_BC.csv"),
                      read_csv("exacerbation_calibration.csv") %>% 
                        mutate(province="CA"))

write_csv(final_result,"master_calibrated_exac.csv")

