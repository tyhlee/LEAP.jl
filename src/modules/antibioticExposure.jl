# include("abstractModule.jl")
# include("../utils")

struct AntibioticExposure <: AntibioticExposure_Module
    hyperparameters::AbstractDict
    parameters::AbstractDict
    AbxOR
end

function process(ag::Agent,anti::AntibioticExposure)
    if !isnothing(anti.parameters[:fixyear])
        if isa(anti.parameters[:fixyear],Number)
            rand(NegativeBinomial(anti.parameters[:θ], antibioticExposure_prob(ag.sex,anti.parameters[:fixyear],anti.parameters)))
        else
            tmp_mu = max(anti.parameters[:midtrends][(ag.cal_year,Int(ag.sex))].rate[1],anti.parameters[:βfloor])
            rand(NegativeBinomial(anti.parameters[:θ], anti.parameters[:θ] / (anti.parameters[:θ] + tmp_mu)))
        end
    else
        rand(NegativeBinomial(anti.parameters[:θ], antibioticExposure_prob(ag.sex,ag.cal_year,anti.parameters)))
    end

end

function process_initial(ag::Agent,anti::AntibioticExposure,cal_born)
    if cal_born < 2001
        rand(NegativeBinomial(anti.parameters[:θ], antibioticExposure_prob(ag.sex,2000,anti.parameters))) 
    else
        if !isnothing(anti.parameters[:fixyear])
            if isa(anti.parameters[:fixyear],Number)
                rand(NegativeBinomial(anti.parameters[:θ], antibioticExposure_prob(ag.sex,anti.parameters[:fixyear],anti.parameters))) 
            else 
                tmp_mu = max( anti.parameters[:midtrends][(cal_born,Int(ag.sex))].rate[1],anti.parameters[:βfloor])
                rand(NegativeBinomial(anti.parameters[:θ], anti.parameters[:θ] / (anti.parameters[:θ] + tmp_mu)))
            end
        else
            rand(NegativeBinomial(anti.parameters[:θ], antibioticExposure_prob(ag.sex,cal_born,anti.parameters))) 
        end
    end
end

# helper function for above
function antibioticExposure_prob(sex::Bool,cal_year,parameters::AbstractDict)
    mu = exp(parameters[:β0] + parameters[:βsex] * sex + parameters[:βcal_year] * cal_year + 
    parameters[:β2005] * (cal_year>2005) +  parameters[:β2005_cal_year] * (cal_year>2005)*cal_year)
    mu = max(mu,parameters[:βfloor]/1000)
    parameters[:θ] / (parameters[:θ] + mu)
end

function random_parameter_initialization!(anti::AntibioticExposure)
    # anti.parameters[:β0] = rand(Normal(anti.hyperparameters[:β0_μ], anti.hyperparameters[:β0_σ]))
    nothing
end
