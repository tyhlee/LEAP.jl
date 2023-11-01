master_birth_estimate = CSV.read(joinpath(dirname(pathof(LEAP)), "processed_data","master_birth_estimate.csv"),DataFrame)
master_population_initial_distribution = CSV.read(joinpath(dirname(pathof(LEAP)), "processed_data","master_initial_pop_distribution_prop.csv"),DataFrame)

master_life_table = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_life_table.csv"),DataFrame);
master_incidence_rate = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_asthma_inc_interpolated.csv"),DataFrame)
master_prevalence_rate = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_asthma_prev_interpolated.csv"),DataFrame)

master_emigration_table = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_emigration_table.csv"),DataFrame)
master_immigration_table = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_immigration_table_modified.csv"),DataFrame)

master_reassessment = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_asthma_assessment.csv"),DataFrame)
master_dx = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_asthma_dx.csv"),DataFrame)
master_mis_dx = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_asthma_mis_dx.csv"),DataFrame)

M3_calibrated_asthma_prev_inc = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_calibrated_asthma_prev_inc_M3.csv"),DataFrame)
exacerbation_calibration  = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","master_calibrated_exac.csv"),DataFrame)

eq5d = CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","eq5d_canada.csv"),DataFrame)
eq5d = groupby(eq5d,[:age,:sex])

abx_mid_trends = groupby(CSV.read(joinpath(dirname(pathof(LEAP)),"processed_data","midtrends.csv"),DataFrame),[:year,:sex])