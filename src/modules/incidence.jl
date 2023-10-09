struct Incidence <: Incidence_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    incidence_table
    prevalence_table
    calibration_table
    min_year
    max_year
    initial_distribution
end

function process(ag::Agent,inc::Incidence,abx)
    tmp_age = min(ag.age,95)
    # max_year = length(inc.incidence_table)
    tmp_year = min(ag.cal_year,inc.max_year)
    if tmp_age < 3
        return false 
    elseif tmp_age < 7
        try 
            rand(Bernoulli(inc.calibration_table[(tmp_year,Int(ag.sex),Int(ag.family_hist),min(ag.num_antibiotic_use,3))][tmp_age-2,"calibrated_inc"]))
        catch 
            println(tmp_year," ",ag.sex, " ", ag.family_hist, " ",ag.num_antibiotic_use, " ", tmp_age)
        end
    else # no effect of Abx beyond 7 years of age
        rand(Bernoulli(inc.calibration_table[(tmp_year,Int(ag.sex),Int(ag.family_hist),0)][tmp_age-2,"calibrated_inc"]))
    end
end

# initialization means prevalence ! ! !
function process_initial(ag::Agent,inc::Incidence)
    tmp_age = min(ag.age,95)-1
    max_year = inc.max_year
    tmp_year = min(ag.cal_year,max_year)-1
    # assume no asthma if age < 3
    if tmp_age < 3
        false
    elseif tmp_age < 7
        rand(Bernoulli(inc.calibration_table[(tmp_year,Int(ag.sex),Int(ag.family_hist),min(ag.num_antibiotic_use,3))][tmp_age-2,"calibrated_prev"]))
    else # no effect of Abx beyond 7 years of age
        rand(Bernoulli(inc.calibration_table[(tmp_year,Int(ag.sex),Int(ag.family_hist),0)][tmp_age-2,"calibrated_prev"]))
    end
end

function process_initial(ag::Agent, inc::Incidence, current_age)    
    # obtain the previous incidence
    min_year = inc.min_year 
    max_year = inc.max_year 
    if current_age==3
        return 3
    else
        find_asthma_age = true
        asthma_age = 3
        tmp_family_hist = Int(ag.family_hist)
        tmp_sex = Int(ag.sex)
        tmp_abx_num = min(ag.num_antibiotic_use,3)
        tmp_year = min(max(ag.cal_year-current_age+asthma_age,min_year),max_year)
        while find_asthma_age && asthma_age < 110
            if rand(Bernoulli(inc.calibration_table[(tmp_year,tmp_sex,tmp_family_hist,(asthma_age < 7 ? min(tmp_abx_num,3) : 0))][asthma_age-2,"calibrated_inc"]))
                return asthma_age
            end
            asthma_age += 1 
            asthma_age = min(asthma_age,95)
            tmp_year += 1
            tmp_year = min(tmp_year,max_year)
        end
        return asthma_age 
    end
end

function prevalence_logit_prob(p::Real,age::Int,param,family_hist,num_abx)
    return p
end

function incidence_logit_prob(p::Real,sex::Bool,age::Int64,cal_year::Int64,ABX::Int64,parameters::AbstractDict,abxlogOR,calibration)
    return p
end


function random_parameter_initialization!(inc::Incidence)
    inc.parameters[:β0] = rand(Normal(inc.hyperparameters[:β0_μ], inc.hyperparameters[:β0_σ]))
end
