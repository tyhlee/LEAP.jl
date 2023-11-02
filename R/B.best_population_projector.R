library(tidyverse)

initial_pop_N <- read_csv("initial_pop_N.csv") %>% 
  filter(year %in% c(2019:2021)) %>% 
  filter(province=="CA") %>% 
  select(-province) %>% 
  rename(REF_DATE = year,
         true_N = N)

pop_projection <- read_csv("~/Dropbox/Asthma_Julia/src/processed_data/pop_all_projected_2020.csv") %>%
  filter(REF_DATE %in% c(2019:2021)) %>% 
  left_join(initial_pop_N,by=c("REF_DATE")) %>% 
  mutate(diff = true_N - all) %>% 
  group_by(Projection_scenario) %>% 
  summarise(n = n(),
            rsq = sqrt(sum((diff)^2))/n) %>% 
  arrange(rsq)


# M3 had the lowest R squred