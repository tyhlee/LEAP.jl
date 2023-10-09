birth_projection = CSV.read(joinpath(dirname(pathof(Asthma_Julia)), "processed_data","birth_projection.csv"),DataFrame)
# change to more intuitive column names
rename!(birth_projection,Dict(:REF_DATE => "calendar_year", :VALUE => "n"))
# correct the data type for n
birth_projection[!,:n] = convert.(Int64,round.(birth_projection.n,digits=0))

master_birth_estimate = CSV.read(joinpath(dirname(pathof(Asthma_Julia)), "processed_data","master_birth_estimate.csv"),DataFrame)
master_population_initial_distribution = CSV.read(joinpath(dirname(pathof(Asthma_Julia)), "processed_data","master_initial_pop_distribution_prop.csv"),DataFrame)

master_life_table = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_life_table.csv"),DataFrame);
master_incidence_rate = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_asthma_inc_interpolated.csv"),DataFrame)
master_prevalence_rate = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_asthma_prev_interpolated.csv"),DataFrame)

emigration_rate = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","pop_emigration_projection_2018.csv"),DataFrame)
emigration_distribution = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","pop_emigration_age_distribution.csv"),DataFrame)
master_emigration_table = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_emigration_table.csv"),DataFrame)

master_immigration_table = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_immigration_table_modified.csv"),DataFrame)

immigration_projection_table = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","pop_immigration_projection_2019.csv"),DataFrame)
immigration_distribution = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","pop_immigration_distribution_2018.csv"),DataFrame)
immigration_rate_per_birth = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","pop_all_projected_2020.csv"),DataFrame)

asthma_initial_distribution = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","asthma_initial_distribution.csv"),DataFrame)
filter!(:sex=> !=("B"),asthma_initial_distribution)
asthma_initial_distribution.sex = (asthma_initial_distribution.sex .== "M")
asthma_initial_distribution.new_prop = rand.(Beta.(asthma_initial_distribution.alpha,asthma_initial_distribution.beta))
asthma_initial_distribution = [filter(:sex => ==(false),asthma_initial_distribution),filter(:sex => ==(true),asthma_initial_distribution)]

master_reassessment = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_asthma_assessment.csv"),DataFrame)
master_dx = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_asthma_dx.csv"),DataFrame)
master_mis_dx = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_asthma_mis_dx.csv"),DataFrame)

M3_calibrated_asthma_prev_inc = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_calibrated_asthma_prev_inc_M3.csv"),DataFrame)
exacerbation_calibration  = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","master_calibrated_exac.csv"),DataFrame)

eq5d = CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","eq5d_canada.csv"),DataFrame)
eq5d = groupby(eq5d,[:age,:sex])

abx_mid_trends = groupby(CSV.read(joinpath(dirname(pathof(Asthma_Julia)),"processed_data","midtrends.csv"),DataFrame),[:year,:sex])