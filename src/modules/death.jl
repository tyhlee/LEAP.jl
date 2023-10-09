# include("abstractModule.jl")
# include("utils")

# prob of death by age and sex
# life_table = CSV.read("../processed_data/life_table.csv",DataFrame);

struct Death <: Death_Module
    parameters
    life_table
end

function process(ag::Agent,d::Death)
    p = d.life_table[ag.cal_year_index][ag.age+1,ag.sex+1]
    if p==1
        return true
    end
    # calibration
    or = p/(1-p)*exp(d.parameters[:β0]+d.parameters[:β1]*ag.cal_year_index+d.parameters[:β2]*ag.age)
    p = max(min(or/(1+or),1),0)
    return rand(Bernoulli(p))
end