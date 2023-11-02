library(tidyverse)

# GOAL study age: 12 - 80

# restrict the EBA cohort to that of the GOAL study
df_goal <- read_rds("private_dataset/df_EBA.rds") %>% 
  filter(age >= 12 & age <=80) %>% 
  filter(visit !=1)

# calculate the annual exacerbation rate
rate <- df_goal$exac %>% mean()
# proportion of asthma control levels:
# controlled, partially controlled, uncontrolled
p_ctls <- c(0.3402965, 0.4737197, 0.1859838)

# calculate the annual exacerbation rate for controlled (c)
# note that rate for partially controlled (pc) is rate_c * 2
# note that rate for uncontrolled (uc) is rate_c * 3
rate_c <- rate/sum(p_ctls*c(1,2,3))
rate_pc <- rate_c*2
rate_uc <- rate_c*3

c(rate_c,rate_pc,rate_uc)
