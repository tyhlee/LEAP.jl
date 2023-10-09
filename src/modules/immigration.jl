# include("abstractModule.jl")
# include("../utils.jl")

# estimated num of newborn + proprotion of male
# birth_project = CSV.read("../processed_data/brith_projection.csv",DataFrame)

struct Immigration <: Immigration_Module
    sex_ratio
    estimate
    age_distribution
    overall_rate
    table
end

function process(sex,age,cal_year,cal_year_index,immi::Immigration)
    Agent(sex,age,cal_year,cal_year_index,true,0,false,nothing,nothing,nothing,[0,0],[zeros(4),zeros(4)],0,false,false)
end
