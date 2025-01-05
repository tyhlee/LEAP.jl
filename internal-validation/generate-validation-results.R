# Import required packages
library(JuliaCall)
library(tidyverse)
library(here)
library(gridExtra)
library(ggpubr)
library(grid)
select <- dplyr::select

# Set up Julia
# Install Julia if necessary: https://julialang.org/downloads/
# Set the path to Julia executable accordingly
julia_path <- "/opt/homebrew/Cellar/julia/1.11.1/bin"
# Initiate a Julia terminal
julia <- JuliaCall::julia_setup(JULIA_HOME = julia_path)
# Import required packages for Julia
julia_command("using JLD2")
julia_command("using JLD")

# Figure settings
fig_setting <- fig_theme <- theme_classic() +
  theme(legend.position = 'top')
# baseline_year = 2001
# chosen_province="CA"
# growth_type = "M3"
# pop_years = c(2001,2010,2015,2020,2025,2030)
# fig_years = c(2001,seq(2005,2030,by=5))
# plot_dims = c(2,3)
# save_fig=T

# Function for generating internal validation results and figures
internal_validator <- function(results,baseline_year=2001,
                               chosen_province='CA',
                               growth_type = "M3",
                               pop_years =  c(2001,2010,2015,2020,2025,2030),
                               fig_years = c(2001,seq(2005,2030,by=5)),
                               plot_dims = c(2,3),
                               save_fig=T,
                               fig_dir = here("internal-validation","figures")){
  
  save_plot <- function(gg_obj,nam,fig.width=10,fig.height=7,fig.dpi=600,save=save_fig){
    if(save){
      ggsave(paste0(fig_dir,'/',nam),gg_obj,width = fig.width,height=fig.height,dpi = fig.dpi)
    } else{
      "plot not saved"
    }
  }
  
  # setting
  max_year <- nrow(results[[1]])
  max_age <- ncol(results[[2]]$immigration[,,1])
  
  longer_df <- function(type){
    tmp_df <- results$outcome_matrix[type][[1]]
    tmp_df_male <- tmp_df[,,2] %>% 
      as.data.frame() %>% 
      mutate(year = row_number()+baseline_year-1) %>% 
      pivot_longer(-year,names_to="age",values_to="n") %>% 
      mutate(age = as.numeric(str_remove(age,"V"))-1) %>% 
      mutate(sex="male")
    tmp_df_female <- tmp_df[,,1] %>% 
      as.data.frame() %>% 
      mutate(year = row_number()+baseline_year-1) %>% 
      pivot_longer(-year,names_to="age",values_to="n") %>% 
      mutate(age = as.numeric(str_remove(age,"V"))-1) %>% 
      mutate(sex="female")  
    
    return(rbind(tmp_df_male,tmp_df_female))
  }
  
  extract_df <- function(type,sex='both',year='all',age='all'){
    if(!is.numeric(sex)){
      sex = ifelse(sex=='male',2,
                   ifelse(sex=='female',1,
                          sex))
    }
    if(year[1]=='all'){
      year <- 1:max_year
    }
    if(age[1]=='all'){
      age <- 1:max_age
    }
    
    if(type %in% c("exacerbation_severity","exacerbation_by_severity")){
      tmp_result <- c()
      for (j in 1:4){
        if (sex=="both"){
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,1,j] + results[[2]][type][[1]][year,age,2,j]
        } else{
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,sex,j]
        }
      }
      
      tmp_result <- tmp_result %>% lapply(.,function(x){
        dim_x <- dim(x)
        matrix(sapply(x,identity),nrow =dim_x[1],dim_x[2])
      })
      
      return(tmp_result)
    } else if(type %in% c("control")){
      tmp_result <- c()
      for (j in 1:3){
        if (sex=="both"){
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,1,j] + results[[2]][type][[1]][year,age,2,j]
        } else{
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,sex,j]
        }
      }
      
      tmp_result <- tmp_result %>% lapply(.,function(x){
        dim_x <- dim(x)
        matrix(sapply(x,identity),nrow =dim_x[1],dim_x[2])
      })
      
      return(tmp_result)
    }
    else{
      if(sex=='both'){
        return(results[[2]][type][[1]][,,1][year,age]+results[[2]][type][[1]][,,2][year,age]) 
      } else{
        return(results[[2]][type][[1]][,,sex][year,age])
      }
    }
    
    
  }
  
  total_n <- sum(results[[2]]$alive)
  
  last_year <- baseline_year+(max_year-1)
  
  #sanity check
  true_N <- results[[1]] %>% rowSums() %>% cumsum()
  alive <- rowSums(extract_df(type = 'alive',sex='both',age='all'))
  death <- rowSums(extract_df(type = 'death',sex='both',age='all'))
  immi <- rowSums(extract_df(type = 'immigration',sex='both',age='all'))
  emi <- rowSums(extract_df(type = 'emigration',sex='both',age='all'))
  

  # 
  # # assumed
  # birth_estimate <- read_csv("src/processed_data/master_birth_estimate.csv") %>% 
  #   filter(year>= baseline_year & year <= baseline_year+nrow(n_generated)-1) %>% 
  #   filter(province==chosen_province) %>% 
  #   filter(projection_scenario %in% c("past",growth_type))
  # true_prop_male <- birth_estimate$prop_male
  # 
  # initial_n <- n_generated$Both[1]
  # 
  # true_prop <- birth_estimate$N/birth_estimate$N[1]

  # Figure 4: death -------------------------------------------------------------------
  # Number of population generated ------------------------------------------
  n_generated <- cbind(Female=extract_df('alive','female',age=1)+extract_df('death','female',age=1)+extract_df("emigration","female",age=1),
                       Male=extract_df('alive','male',age=1)+extract_df('death','male',age=1)+extract_df("emigration","male",age=1)) %>%
    as.data.frame() %>%
    mutate(year=row_number()+baseline_year-1,
           Both=Female+Male)
  
  life_table <- read_csv(here("src","processed_data","master_life_table.csv")) %>%
    filter(year>=baseline_year & year <= last_year & province==chosen_province) %>% 
    # rename(Female=prob_death_female,
    #        Male=prob_death_male,
    #        Age = age) %>% 
    #   pivot_longer(cols=-Age,values_to='prop',names_to='sex') %>% 
    mutate(Type='True') %>% 
    mutate(sex = ifelse(sex=="M","Male","Female")) %>% 
    select(year,age,sex,prob_death,Type) %>% 
    rename(prop=prob_death)
  
  death_female <- extract_df(type = 'death',sex='female',age='all')/cbind(n_generated$Female,extract_df(type = 'alive',sex='female',age='all')[,-max_age])
  death_male <- extract_df(type = 'death',sex='male',age='all')/cbind(n_generated$Male,extract_df(type = 'alive',sex='male',age='all')[,-max_age])
  
  # death_female <- extract_df(type = 'death',sex='female',age='all')/extract_df(type = 'alive',sex='female',age='all')
  # death_male <- extract_df(type = 'death',sex='female',age='all')/extract_df(type = 'alive',sex='male',age='all')
  
  tmp_death <- c()
  for(i in 1:(max_year)){
    tmp_death[[i]] <- data.frame(year=baseline_year+i-1,
                                 age=1:max_age,
                                 Female=death_female[i,],
                                 Male=death_male[i,])
  }
  df_death <- do.call(rbind,tmp_death) %>% 
    mutate(age=age-1) %>% 
    pivot_longer(cols=-c(1:2),values_to='prop',names_to='sex') %>% 
    mutate(Type='Observed')
  
  death_plotter <- function(chosen_year){
    ggplot(data=rbind(df_death,life_table) %>% 
             filter(age<=80) %>% 
             filter(year==chosen_year),
           aes(x=age,y=prop,color=Type,linetype=Type)) +
      geom_line(linewidth=2)+
      facet_grid(.~sex)+
      # scale_linetype_manual(values=c("twodash", "dotted"))+
      ylab("")+
      xlab('') +
      ggtitle(chosen_year)+
      fig_theme +
      theme(legend.title=element_blank())+
      scale_color_manual(values=c("grey","black"),c("Observed","True"))
  }
  # pop_years <- c(2000,2004,2008,2012,2015,2018)
  
  fig.death <- do.call("ggarrange", c(lapply(pop_years,function(x){death_plotter(x)}),common.legend=T,legend='none')) %>% 
    annotate_figure(left = textGrob("Probabilty of dying between ages x and x+1", rot = 90, vjust = 1,
                                    gp=gpar(fontsize=15)),
                    bottom = textGrob("Age (year)", gp = gpar(fontsize = 15)))
  save_plot(fig.death,'fig_4.jpeg')
  
  
  # Figure 5: Population pyramid -------------------------------------------------------
  n_alive_female <- extract_df('alive','female',age='all') +
    extract_df('death','female',age='all') + extract_df("emigration","female",age='all')
  n_alive_male <- extract_df('alive','male',age='all')   +
    extract_df('death','male',age='all') + extract_df("emigration","male",age='all')
  n_alive_both <- n_alive_female + n_alive_male
  
  pop_est <- read_csv("src/processed_data/initial_pop_distribution.csv") %>%
    filter(province==chosen_province) %>%
    filter(year >= baseline_year) %>%
    mutate(sex=ifelse(sex,"Male","Female"))

  pop_CA_BC <- read_rds(here("src","processed_data","pop_projection_BC_CA.rds")) %>%
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
    filter(year <= last_year)

  comparsion <- function(chosen_year=baseline_year){
    chosen_year_index <- chosen_year - baseline_year + 1
    est <- n_alive_both[chosen_year_index,]
    # est[101] <- sum(est[101:length(est)])
    #0.1206846
    true <- pop %>% 
      filter(year==chosen_year) %>% 
      group_by(year,age) %>% 
      summarise(n=sum(n)) %>%
      ungroup() 
    # filter(age<=100)
    
    if(chosen_year==2000){
      true <- true %>% 
        rbind(data.frame(year=2000,age=100,n=0))
    } 
    
    ratio <- true$n[1]/est[1]
    true <- cbind(true,estimated=est[1:101]*ratio) %>% 
      filter(age<100)
    ggplot(data=true %>% 
             rename(True=n,
                    Observed=estimated) %>% 
             pivot_longer(cols=-c(year:age),names_to='Type',values_to="prop"),
           aes(x=age,y=prop,linetype=Type,col=Type))+
      geom_line(linewidth=2) +
      ggtitle(chosen_year)+
      xlab("")+
      ylab("") +
      fig_theme +
      theme(legend.title=element_blank(),
            legend.position ='none') +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values=c("grey","black"),c("Observed","True"))
  }
  
  fig.pryamid <- do.call("ggarrange", c(lapply(pop_years,function(x){comparsion(x)}), 
                                        common.legend=T,legend='none')) %>% 
    annotate_figure(left = textGrob("Number of population", rot = 90, vjust = 1,
                                    gp=gpar(fontsize=15)),
                    bottom = textGrob("Age (year)", gp = gpar(fontsize = 15)))
  save_plot(fig.pryamid,'fig_5.jpeg')
  
  # Family history ----------------------------------------------------------
  rate_fam_history <-  sum(extract_df("family_history","both"))/sum(results[[1]])
  sd <- 1.96*sqrt(sum(extract_df("family_history","both"))/sum(results[[1]])^2)
  rate_fam_history+sd
  rate_fam_history-sd
  
  # should be around 0.2927242
  if(abs(rate_fam_history-0.2927242) < 0.01 ){
    fam_result <- "family history is validated"
  } else{
    fam_result <- "family history is NOT validated internally"
  }
  
  print(fam_result)
  
  # Figure 6: Antibiotic exposure ------------------------------------------------------
  # count_model <- read_rds("data-raw/count_model_BC.rds")
  # count_data <- read_rds("data-raw/Abx_count_data_BC.rds") %>%
  #   mutate(rate=N_Abx/N*1000) %>%
  #   select(year,sex,rate) %>%
  #   mutate(sex=ifelse(sex==0,"Female","Male"))
  count_data <- read_csv("src/processed_data/InfantAbxBC.csv")
  
  df_ABE <- rbind(data.frame(Year=1:max_year,
                             n=extract_df(type = 'antibiotic_exposure',sex='male',age=1),
                             N = extract_df(type = 'alive',sex='male',age=1)+
                               extract_df(type = 'death',sex='male',age=1)+
                               extract_df(type = 'emigration',sex='male',age=1)) %>% 
                    mutate(sex=1),
                  data.frame(Year=1:max_year,
                             n=extract_df(type = 'antibiotic_exposure',sex='female',age=1),
                             N = extract_df(type = 'alive',sex='female',age=1)+
                               extract_df(type = 'death',sex='female',age=1)+
                               extract_df(type = 'emigration',sex='female',age=1)) %>% 
                    mutate(sex=0))%>% 
    mutate(rate=n/N*1000) %>% 
    mutate(sex=ifelse(sex==0,"Female","Male")) %>% 
    mutate(Year=Year+baseline_year-1)
  
  # scaler <- round(max(1/(df_ABE$rate/df_ABE$n)))
  
  ggplot(data=df_ABE) +
    geom_line(aes(y=rate,x=as.numeric(Year),col=sex))+
    geom_line(data=count_data,aes(y=rate,x=as.numeric(year),col=sex),linetype='dashed')+
    # geom_col(aes(y=n/scaler,x=as.numeric(Year),fill="#bdbdbd")) +
    # scale_y_continuous(sec.axis=sec_axis(~(.)*scaler,name="Number of antibotics")) +
    theme_classic() +
    # scale_linetype_identity(name=NULL,labels=c("Rate"),guide='legend')+
    # scale_fill_identity(name=NULL,labels="Number",guide='legend')+
    ylab("Number of courses of antibiotics in the first year of life (per 1,000)")+
    scale_y_continuous(breaks=seq(0,1000,by=100)) +
    scale_x_continuous(breaks=fig_years) +
    geom_hline(yintercept = 50,col='purple',linetype='dashed')+
    xlab('Year') +
    fig_setting +
    theme(legend.title = element_blank())-> fig.ABE
  
  save_plot(fig.ABE,'fig_6.jpeg')
  
  # Figure 7: Asthma prevalence -------------------------------------------------------
  
  simulation_max_year <- baseline_year+max_year-1
  
  asthma_inc <- read_csv("src/processed_data/master_asthma_prev_inc.csv") %>% 
    mutate(province=chosen_province) %>% 
    select(year,age,province,sex,inc) %>% 
    rename(value = inc) %>% 
    mutate(sex = ifelse(sex==0,"F","M")) %>% 
    mutate(observed=NA,
           n=NA,
           N=NA) %>% 
    filter(year <= simulation_max_year) %>% 
    filter( year >= baseline_year)
  
  asthma_prev <- read_csv("src/processed_data/master_asthma_prev_inc.csv") %>% 
    mutate(province=chosen_province) %>% 
    select(year,age,province,sex,prev) %>% 
    rename(value = prev) %>% 
    mutate(sex = ifelse(sex==0,"F","M")) %>% 
    mutate(observed=NA,
           n=NA,
           N=NA) %>% 
    filter(year <= simulation_max_year) %>% 
    filter( year >= baseline_year)
  
  for(i in 1:nrow(asthma_inc)){
    tmp <- asthma_inc[i,]
    tmp_age <- tmp$age+1
    tmp_sex <- as.numeric(tmp$sex=="M")+1
    tmp_year <- tmp$year-baseline_year +1
    asthma_inc$n[i] <- extract_df(type = 'asthma_incidence',sex=tmp_sex,year=tmp_year,age=tmp_age)
    asthma_prev$n[i] <- extract_df(type = 'asthma_prevalence',sex=tmp_sex,year=tmp_year,age=tmp_age)
    asthma_inc$N[i] <- asthma_prev$N[i] <- sum(extract_df(type = 'alive',sex=tmp_sex,age=tmp_age,year=tmp_year) +
                                                 extract_df(type = 'death',sex=tmp_sex,age=tmp_age,year=tmp_year) +
                                                 extract_df(type = 'emigration',sex=tmp_sex,age=tmp_age,year=tmp_year))
  }
  
  asthma_prev <- asthma_prev %>% 
    mutate(observed=n/N) %>% 
    rename(true=value) %>% 
    pivot_longer(5:6,names_to="Type") %>% 
    mutate(CI_lower = ifelse(Type=='true',NA,value-1.96*sqrt(value*(1-value)/N)),
           CI_upper = ifelse(Type=='true',NA,value+1.96*sqrt(value*(1-value)/N))) %>% 
    mutate(value = value*1000,
           CI_lower = CI_lower*1000,
           CI_upper = CI_upper * 1000)
  
  
  ggplot(asthma_prev %>% 
           filter(!is.na(CI_lower)) %>% 
           filter(year %in% c(2020,2025,2030,2035,2040,2050)) %>% 
           mutate( year = as.factor(year)),
         aes(y=value,x=age,colour=year)) +
    geom_line() +
    xlim(c(0,60))+
    facet_grid(.~sex)
  
  asthma_inc <- asthma_inc %>% 
    mutate(observed=n/N)%>% 
    rename(true=value) %>% 
    pivot_longer(5:6,names_to="Type") %>% 
    mutate(CI_lower = ifelse(Type=='true',NA,value-1.96*sqrt(value*(1-value)/N)),
           CI_upper = ifelse(Type=='true',NA,value+1.96*sqrt(value*(1-value)/N))) %>% 
    mutate(value = value*1000,
           CI_lower = CI_lower*1000,
           CI_upper = CI_upper * 1000)
  
  
  
  asthma_plotter <- function(chosen_year,inc=T){
    
    if(inc){
      tmp_df <- asthma_inc
      
      ggplot(data=tmp_df %>% 
               filter(year == chosen_year) %>% 
               mutate(sex = ifelse(sex=="F","Female","Male"))) +
        geom_line(aes(y=value,x=age,col=Type,linetype=Type),linewidth=2)+
        facet_grid(.~sex) + 
        theme_classic() +
        fig_theme+
        ylab("")+
        ylim(c(0,100))+
        xlim(c(0,65))+
        theme(legend.position='none')+
        ggtitle(chosen_year)+
        xlab('') +
        scale_color_manual(values=c("grey","black"),c("Observed","True"))+
        theme(text = element_text(size=15))
    } else{
      tmp_df <- asthma_prev
      
      ggplot(data=tmp_df %>% 
               filter(year == chosen_year) %>% 
               mutate(sex = ifelse(sex=="F","Female","Male"))) +
        geom_line(aes(y=value,x=age,col=Type,linetype=Type),linewidth=2)+
        facet_grid(.~sex) + 
        theme_classic() +
        fig_theme+
        ylab("")+
        ylim(c(0,250))+
        xlim(c(0,65))+
        theme(legend.position='none')+
        ggtitle(chosen_year)+
        xlab('') +
        scale_color_manual(values=c("grey","black"),c("Observed","True")) +
        theme(text = element_text(size=15))
    }
  }
  
  fig.asthma.prev <- do.call("ggarrange", c(lapply(pop_years,function(x){asthma_plotter(x,F)}), 
                                            common.legend=T,legend='none')) %>% 
    annotate_figure(left = textGrob("Asthma prevalence (per 1,000)", rot = 90, vjust = 1,
                                    gp=gpar(fontsize=20)),
                    bottom = textGrob("Age (year)", gp = gpar(fontsize = 20)))
  save_plot(fig.asthma.prev,'fig_7.jpeg')
  
  # Asthma prevalence OR ---------------------------------------------------------------
  options(dplyr.summarise.inform = FALSE)
  
  OR_abx_calculator <- function(age,dose,params = c(1.711+0.115,-0.225,0.053)){
    ifelse(dose==0,1,exp(sum(params*c(1,pmin(age,7),pmin(dose,3)))))
  }
  
  OR_fam_calculator <- function(age,fam_hist,params=c(log(1.13),(log(1.13)+log(2.4))/2-log(1.13))){
    ifelse(age<3 | fam_hist==0 | age > 7,1,exp(params[1] + params[2]*(pmin(age,5)-3)))
  }
  
  OR_risk_factor_calculator <- function(fam_hist,age,dose,params=list(c(log(1.13),(log(1.13)+log(2.4))/2-log(1.13)),c(1.711+0.115,-0.225,0.053))){
    ifelse(age<3,1,exp(log(OR_fam_calculator(age,fam_hist,params[[1]])) + log(OR_abx_calculator(age,dose,params[[2]]))))
  }
  
  true_inc_table <- true_prev_table <- rbind(expand.grid(fam_history = c(0,1),
                                                         abx_exposure = c(0:3),
                                                         year = c(2001:2030),
                                                         age = 3:7,
                                                         sex=c(0,1),
                                                         true_OR = NA),
                                             expand.grid(fam_history = c(0,1),
                                                         abx_exposure = 0,
                                                         year = c(2001:2030),
                                                         age = 8:65,
                                                         sex=c(0,1),
                                                         true_OR = NA)) %>% 
    as.data.frame()
  
  true_prev_table$true_OR <- apply(true_prev_table,1,FUN=function(x){
    OR_risk_factor_calculator(x[1],x[4],x[2],params=list(c(log(1.13),(log(1.13)+log(2.4))/2-log(1.13)),c(1.711+0.115,-0.225,0.053)))
  })
  
  true_inc_table$true_OR <- apply(true_inc_table,1,FUN=function(x){
    OR_risk_factor_calculator(x[1],x[4],x[2],params=list(c(log(1.13),0.3619942),c(1.711+0.115,-0.2920745,0.053)))
  })
  
  raw_prev_table <- data.frame(results$outcome_matrix$asthma_prevalence_contingency_table)
  colnames(raw_prev_table) <- c("year","sex","fam_history","abx_exposure","age","n_asthma","n_no_asthma")
  raw_prev_table %>% 
    filter(age>2) %>% 
    filter(age <= 65) %>% 
    mutate(prevalence = n_asthma/ (n_asthma + n_no_asthma)) %>% 
    left_join(true_prev_table,by=c("year","sex","age","fam_history","abx_exposure")) %>% 
    arrange(year,sex,age,fam_history,abx_exposure) %>% 
    group_by(year,sex,age) %>% 
    group_split() %>% 
    lapply(.,function(x){
      ref_row <- x[1,]
      tmp_nrow <- nrow(x)
      tmp_age <- ref_row$age
      tmp_OR <- tmp_OR_CI_u <- tmp_OR_CI_l <- c()
      tmp_OR[[1]] <- tmp_OR_CI_u[[1]] <- tmp_OR_CI_l[[1]] <- 1
      
      if(tmp_age < 8){
        for(i in 1:(tmp_nrow-1)){
          if(any(is.na(c(ref_row$n_no_asthma,ref_row$n_no_asthma, x$n_asthma[i+1],x$n_no_asthma[i+1])))) {
            tmp_OR[[i+1]] <- tmp_OR_CI_u[[i+1]] <- tmp_OR_CI_l[[i+1]] <-  NA
          } else{
            calc_OR <- fisher.test(matrix(c(ref_row$n_no_asthma,ref_row$n_asthma,x$n_no_asthma[i+1], x$n_asthma[i+1]),ncol = 2,byrow = T))
            tmp_OR[[i+1]] <- calc_OR$estimate
            tmp_OR_CI_l[[i+1]] <- calc_OR$conf.int[1]
            tmp_OR_CI_u[[i+1]] <- calc_OR$conf.int[2]
          }
        }
        x$observed_OR <- unlist(tmp_OR)
        x$observed_OR_CI_u <- unlist(tmp_OR_CI_u)
        x$observed_OR_CI_l <- unlist(tmp_OR_CI_l)
        
        return(x)
      } else {
        x %>% 
          mutate(abx_exposure=0) %>% 
          group_by(year,sex,fam_history,abx_exposure,age) %>% 
          summarise(n_asthma = sum(n_asthma),
                    n_no_asthma = sum(n_no_asthma),
                    # prob = sum(prob,na.rm=T),
                    # true_prev = sum(true_prev,na.rm=T),
                    true_OR = sum(true_OR,na.rm=T)) %>%
          mutate(prevalence = n_asthma/n_no_asthma) %>% 
          select(year:n_no_asthma,prevalence,
                 # prob,true_prev,
                 true_OR) -> tmp_x
        
        if(any(is.na(c(tmp_x$n_asthma,tmp_x$n_no_asthma)))){
          tmp_x$observed_OR <- c(1,NA)
          tmp_x$observed_OR_CI_u <- c(1,NA)
          tmp_x$observed_OR_CI_l <- c(1,NA)
          return(tmp_x)
        } else{
          calc_OR <- fisher.test(matrix(c(tmp_x$n_no_asthma[1],tmp_x$n_asthma[1], 
                                          tmp_x$n_no_asthma[2],tmp_x$n_asthma[2]),ncol = 2,byrow = T))
          tmp_x$observed_OR <- c(1,calc_OR$estimate)
          tmp_x$observed_OR_CI_u <- c(1,calc_OR$conf.int[2])
          tmp_x$observed_OR_CI_l <- c(1,calc_OR$conf.int[1])
          
          return(tmp_x %>% ungroup())
        }
      }
    }) %>% 
    do.call(rbind,.)  %>% 
    group_by(year,sex,age) %>% 
    mutate(total_n = n_asthma+n_no_asthma,
           group_n = sum(total_n),
           observed_prob = total_n/group_n) %>% 
    ungroup() %>% 
    mutate(exp_abs_diff_log_OR = exp(abs(log(true_OR) - log(observed_OR))),
           # prev_checker = abs(prevalence-true_prev)/true_prev*100,
           within_CI = true_OR <= observed_OR_CI_u & true_OR >= observed_OR_CI_l) -> prev_table
  
  prev_OR <- prev_table %>%
    mutate(key = paste0(year,"_",age,"_",sex)) %>%
    filter(true_OR!=1) %>% 
    group_by(year,age,sex) %>%
    mutate(ratio_OR = log(true_OR/observed_OR),
           exp_abs_diff_log_OR = exp(abs(log(true_OR) - log(observed_OR))),
           within_CI = true_OR <= observed_OR_CI_u & true_OR >= observed_OR_CI_l) %>%
    filter(!(is.na(exp_abs_diff_log_OR) | is.infinite(exp_abs_diff_log_OR)))
  
  summary(prev_OR$ratio_OR)
  
  quantile(prev_OR$ratio_OR,c(0.01,0.99))
  
  prev_OR %>% 
    filter(ratio_OR < -0.40 | ratio_OR > 0.70) %>% 
    View()
  
  # Figure 8: Control -----------------------------------------------------------------
  n_control <- extract_df("control",'both') %>% 
    lapply(.,function(x){
      rowSums(x)
    })
  
  df_control <- data.frame(Year=1:max_year,
                           n1=n_control[[1]],
                           n2=n_control[[2]],
                           n3=n_control[[3]]) %>% 
    mutate(N=n1+n2+n3,
           rate1=n1/N,
           rate2=n2/N,
           rate3=n3/N)
  
  ggplot(data=df_control %>% 
           select(Year,rate1:rate3) %>% 
           mutate(Year = Year+baseline_year) %>% 
           pivot_longer(cols=-Year,names_to="Type",values_to='prop') %>% 
           mutate(Type = case_when(Type=='rate1' ~ 'well-controlled',
                                   Type=='rate2' ~ 'partially controlled',
                                   TRUE ~ 'uncontrolled'),
                  Type = factor(Type,c("well-controlled","partially controlled","uncontrolled"),
                                c("well-controlled","partially controlled","uncontrolled"))),
         aes(x=Year,y=prop,color=Type)) +
    geom_line(linewidth=2,alpha=0.7)+
    fig_theme+
    ylim(0,0.50)+
    geom_hline(yintercept = 0.47,col='green',linetype='dashed',linewidth=2)+
    geom_hline(yintercept = 0.35,col='red',linetype='dashed',linewidth=2)+
    geom_hline(yintercept = 0.18,col='blue',linetype='dashed',linewidth=2)+
    # (0.3476309, 0.4748130, 0.1775561)
    ylab("Proportion of time spent in each control level")+
    xlab('Year') +
    scale_x_continuous(breaks=fig_years)+
    fig_setting +
    theme(text=element_text(size=20))+
    theme(legend.position = 'none')-> fig.control
  
  save_plot(fig.control,'fig_8.jpeg')

  # Figure 9: Exacerbation -------------------------------------------------
  df_exac <- data.frame(Year=1:max_year,
                        n=rowSums(extract_df(type = 'exacerbation',sex='both',age='all')),
                        N = rowSums(extract_df(type = 'asthma_prevalence',sex='both',age='all'))) %>%
    mutate(rate=n/N)
  
  n_severity <- extract_df("exacerbation_by_severity",'both') %>% 
    lapply(.,function(x){
      rowSums(x)
    })
  
  df_severity <- data.frame(Year=1:max_year,
                            n1=n_severity[[1]],
                            n2=n_severity[[2]],
                            n3=n_severity[[3]],
                            n4=n_severity[[4]]) %>% 
    mutate(N=n1+n2+n3+n4,
           rate1=n1/N,
           rate2=n2/N,
           rate3=n3/N,
           rate4=n4/N)
  
  ggplot(data=df_severity %>% 
           select(Year,rate1:rate4) %>% 
           mutate(Year = Year+baseline_year) %>% 
           pivot_longer(cols=-Year,names_to="Type",values_to='prop') %>% 
           mutate(Type = case_when(Type=='rate1' ~ 'mild',
                                   Type=='rate2' ~ 'moderate',
                                   Type == 'rate3' ~ "severe",
                                   TRUE ~ 'very severe')),
         aes(x=Year,y=prop,color=Type)) +
    geom_line(linewidth=2,alpha=0.7)+
    fig_theme+
    theme(legend.title=element_blank())+
    ylim(0,0.5)+
    scale_color_discrete(c("rate1","rate2","rate3","rate4"),
                         c("red","green","blue","purple"))+
    geom_hline(yintercept = 0.495,col='red',linetype='dashed',linewidth=2)+
    geom_hline(yintercept = 0.195,col='green',linetype='dashed',linewidth=2)+
    geom_hline(yintercept = 0.283,col='blue',linetype='dashed',linewidth=2)+
    geom_hline(yintercept = 0.026,col='purple',linetype='dashed',linewidth=2)+
    ylab("Proportion of exacerbation severity")+
    xlab('Year') +
    scale_x_continuous(breaks=fig_years)+
    fig_setting +
    theme(text=element_text(size=20))+
    theme(legend.title=element_blank())-> fig.exac_sev
  
  save_plot(fig.exac_sev,'fig_9.jpeg')
  
  
  # Figure 10: Hospitalization -------------------------------------------
  df_cihi <- read_rds(paste0("R/public_dataset/asthma_hosp/",
                             chosen_province,"/tab1.rds"))$rate %>% 
    filter(fiscal_year >= min(2001,baseline_year)) %>% 
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
    select(-type)
  
  
  tmp_cihi <- df_cihi%>% filter(year == max(df_cihi$year))
  impute_years <- (max(df_cihi$year)+1):(max_year+baseline_year)
  
  for(i in impute_years){
    df_cihi <- rbind(df_cihi,
                     tmp_cihi %>% 
                       mutate(year=i))
  }
  
  
  df_alive <- longer_df("alive")
  df_death <- longer_df("death")
  df_emi <- longer_df("emigration")
  df_asthma <- longer_df("asthma_prevalence") %>% 
    rename(n_asthma = n)
  df_pop <- df_alive %>% 
    left_join(df_death,by=c("year","age","sex")) %>% 
    left_join(df_emi,by=c("year","age","sex")) %>% 
    mutate(total_n = n+n.x+n.y) %>% 
    select(year,age,sex,total_n) %>% 
    left_join(df_asthma,by=c("year","age","sex"))
  
  df_hosp <- longer_df("exacerbation_hospital") %>% 
    left_join(df_pop,by=c("year","age","sex")) %>% 
    mutate(rate = n/total_n*100000) %>% 
    filter(age>=3) %>% 
    mutate(sex=as.numeric(sex=="male")) %>% 
    mutate(age= ifelse(age>90,90,age)) %>% 
    group_by(year,sex,age) %>% 
    summarise(n=sum(n),
              total_n = sum(total_n),
              n_asthma = sum(n_asthma)) %>% 
    mutate(rate = n/total_n*100000,
           rate_per_asthma = n/n_asthma*100000) %>% 
    left_join(df_cihi,by=c("year","age","sex")) %>% 
    # filter(year<=2017) %>% 
    mutate(rate_upper = rate+qnorm(0.975)*sqrt(n/total_n^2)*100000,
           rate_lower =rate-qnorm(0.975)*sqrt(n/total_n^2)*100000,
           within_CI = true_rate <= rate_upper & true_rate >= rate_lower) %>% 
    mutate(ratio = rate/true_rate) %>% 
    filter(age <= 65)
  
  df_cihi_count <- read_rds(paste0("R/public_dataset/asthma_hosp/",
                                   chosen_province,"/tab1.rds"))$count %>% 
    filter(fiscal_year >= min(2001,baseline_year)) %>% 
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
    filter(!(age %in% c(0:2))) %>% 
    group_by(year,sex) %>% 
    summarise(count=sum(true_rate))
  
  df_cihi_N <- read_rds(paste0("R/public_dataset/asthma_hosp/",
                               chosen_province,"/tab1.rds"))$N %>% 
    filter(fiscal_year >= min(2001,baseline_year)) %>% 
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
    filter(!(age %in% c(0:2))) %>% 
    group_by(year,sex) %>% 
    summarise(N=sum(true_rate))
  
  df_cihi_rate <- df_cihi_count %>% 
    left_join(df_cihi_N,by=c("year","sex")) %>% 
    mutate(true_rate = count/N*100000)
  
  overall_hosp <- df_hosp %>% 
    group_by(year,sex) %>% 
    summarise(n=sum(n),
              total_n= sum(total_n)) %>% 
    mutate(rate=n/total_n*100000) %>% 
    left_join(df_cihi_rate %>% 
                select(year,sex,true_rate),
              by=c("year","sex")) %>% 
    mutate(rate_upper = rate+qnorm(0.975)*sqrt(n/total_n^2)*100000,
           rate_lower =rate-qnorm(0.975)*sqrt(n/total_n^2)*100000) %>% 
    mutate(sex = ifelse(sex==1,"Males","Females")) %>% 
    mutate(ratio = true_rate/rate,
           log_ratio = log(ratio))
  
  hosp_summary <- (overall_hosp %>% filter(year<=2017))$log_ratio %>% summary()
  
  
  ggplot(overall_hosp %>% 
           select(1:2,5:6) %>% 
           pivot_longer(cols=-c(1,2),names_to="Type",values_to="rate") %>% 
           mutate(type=ifelse(Type=='rate',"Observed","True")),
         aes(x=year,y=rate,col=Type,linetype=Type)) +
    geom_line(linewidth=2) +
    theme_classic() +
    facet_grid(.~sex)+
    ylim(0,120) +
    theme(legend.position='none') +
    scale_color_manual(values=c("grey","black"),c("Observed","True")) +
    ylab("Rate of very severe exacerabtions per 100,000") +
    theme(text=element_text(size=20))+
    scale_x_continuous(breaks=fig_years)+
    xlab("Year") -> fig.asthma.hosp
  
  save_plot(fig.asthma.hosp,'fig_10.jpeg')
  
  return(list(fig.death,
              fig.pryamid,
              fam_result,
              fig.ABE,
              fig.asthma.prev,
              prev_OR,
              fig.control,
              fig.exac_sev,
              fig.asthma.hosp,
              hosp_summary))
}

# Load simulation results
results <- julia_eval('JLD2.load("internal-validation/LEAP_output.jld")["output"]') 

# Generate internal validation results
IV_results <- internal_validator(results)

prev_OR <- IV_results[[6]]
summary(prev_OR$ratio_OR)
quantile(prev_OR$ratio_OR,c(0.01,0.99))
# Look at extreme cases:
prev_OR %>% 
  filter(ratio_OR < -0.40 | ratio_OR > 0.70) %>% 
  View()
