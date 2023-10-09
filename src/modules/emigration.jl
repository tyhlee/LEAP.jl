# include("abstractModule.jl")
# include("../utils.jl")

# estimated num of newborn + proprotion of male
# birth_project = CSV.read("../processed_data/brith_projection.csv",DataFrame)

struct Emigration <: Emigration_Module
    projected_rate
    age_distribution
    table
end

# function process(ag_age::Int,b::Emigration)
#     rand(Bernoulli(b.projected_rate*b.age_distribution.percentage[searchsortedfirst(b.age_distribution.age_upper,ag_age)]))
# end

function process(cal_year_index::Int,age::Int,sex::Bool,b::Emigration)
    if age == 0
        return false
    else 
        return rand(Bernoulli(b.table[cal_year_index][min(age,100),sex+3]))
    end
end