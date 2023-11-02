library(tidyverse)
library(here)

death_cause <- read_csv('public_data/13100782.csv') %>% 
  select(1:2,4:6,VALUE) 
colnames(death_cause) <- gsub(' ','_',colnames(death_cause))
colnames(death_cause)[5] <- 'cause_death_ICD10'
look <- death_cause$cause_death_ICD10 %>% unique()
asthmas <-look[str_detect(look,"asthma|Asthma")]

death_cause <- death_cause %>% 
  mutate(Sex = str_sub(Sex,1,1),
         Age_group = if_else(Age_group=='Total, all ages','all',
                             if_else(Age_group=='Under 1 year','0-1',
                                     gsub(' ','',gsub('years','',gsub('to','-', Age_group))))))

death_cause_asthma <- death_cause %>% 
  filter(cause_death_ICD10 %in% c("Asthma [J45]","Status asthmaticus [J46]") &
           Age_group != 'Age, not stated') %>% 
  group_by(REF_DATE,GEO,Age_group,Sex) %>% 
  summarise(value = sum(VALUE)) %>% 
  ungroup()


death_cause_asthma %>% 
  filter(Age_group=='all' & Sex=='B') %>%
  left_join(death_cause %>% 
              filter(cause_death_ICD10 == "Total, all causes of death [A00-Y89]") %>%  
              select(-cause_death_ICD10),
            by=c("REF_DATE","GEO","Age_group","Sex")) %>% 
  mutate(prop_asthma_death = value/VALUE) -> death_due_to_asthma


summary(death_due_to_asthma$prop_asthma_death)
# max of annual proportions of asthma-related deaths is 0.14% 
# ignore

# comparison to COPD
COPD <-look[str_detect(look,"obstructive pulmonary")]

death_cause_COPD <- death_cause %>% 
  filter(cause_death_ICD10 %in% "Other chronic obstructive pulmonary disease [J44]" &
           Age_group != 'Age, not stated') %>% 
  group_by(REF_DATE,GEO,Age_group,Sex) %>% 
  summarise(value = sum(VALUE)) %>% 
  ungroup()


death_cause_COPD %>% 
  filter(Age_group=='all' & Sex=='B') %>% 
  left_join(death_cause %>% 
              filter(cause_death_ICD10 == "Total, all causes of death [A00-Y89]") %>%  
              select(-cause_death_ICD10),
            by=c("REF_DATE","GEO","Age_group","Sex")) %>% 
  mutate(prop_death = value/VALUE) -> death_due_to_COPD

# conclusion: no need to adjust for deaths due to asthma