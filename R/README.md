| Number | File name | Description | Output |
|-|:-:|:------:|:---:|
| 1 | birth.R | estimated and projected number of births and the proportion of males (vs. females) by year, age, province, and projection scenario | master_birth.csv |
| 2 | initial_population.R | estimated and projected number of people, proportion of people relative to infants (aged less than 1 year), and proportion of males (vs. females) by year, age, province, and projection scenario | master_initial_pop_distribution_prop.csv |
| 3 | mortality.R | estimated life table by year, age, sex, and province |life_table.csv |
| 4 | calibrate_lifetable_immigration_emigration.R | Projected life table, net immigrants and emigrants, by year, sex, age  | master_life_table.csv,master_emigration_table.csv, master_immigration_table.csv, master_immigration_table_modified.csv |
| 5 |asthma_occurrence_crude.R | Interpolated asthma incidence and prevalence by year, sex, age, and province | master_asthma_inc_interpolated.csv, master_asthma_prev_interpolated.csv|
| 6 | asthma_dx.R | Diagnosis and reassessment of asthma by year, sex, age, province | master_asthma_assessment.csv, master_asthma_dx.csv, master_asthma_mis_dx.csv |
| 7 | asthma_inc_prev_calibration.R, calibration_helper_function.R| Calibrate the effect of risk factors on asthma incidence and prevalence | master_calibrated_asthma_prev_inc_M3.csv |
| 8 | exacerbation_calibration.R | Calibrate the asthma hospitalization rate | master_calibrated_exac.csv|
| 9 | asthma_control.pdf | Analysis of asthma control based on the Economic Burden of Asthma data  | asthma control model |
| 10 | utility.R | Generate the utility values for the general population by sex and age | eq5d_canada.csv |
| A | asthma_mortality_adjustment.R  | Investigate whether mortality adjustment is necessary for asthma patients | |
| B | best_population_project.R| Determine the best population projection scenario | | 
