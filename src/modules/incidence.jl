struct Incidence <: Incidence_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    parameters_prev::Union{AbstractDict,Nothing}
    min_year
    max_year
    max_age
end

function process(ag::Agent,inc::Incidence)
    tmp_age = min(ag.age,inc.max_age)
    if tmp_age < 3
        return false 
    elseif tmp_age ==3
        return rand(Bernoulli(prevalence_equation(ag.sex,tmp_age,ag.cal_year,ag.family_hist,ag.num_antibiotic_use,inc.parameters_prev,inc.max_year)))
    else
        return rand(Bernoulli(incidence_equation(ag.sex,tmp_age,ag.cal_year,ag.family_hist,ag.num_antibiotic_use,inc.parameters,inc.max_year)))
    end
end

# initialization means prevalence ! ! !
function process_initial(ag::Agent,inc::Incidence)
    tmp_age = min(ag.age-1,inc.max_age)
    tmp_year = ag.cal_year-1
    # assume no asthma if age < 3
    if tmp_age < 3
        return false
    else # no effect of Abx beyond 7 years of age
        return rand(Bernoulli(prevalence_equation(ag.sex,tmp_age,tmp_year,ag.family_hist,ag.num_antibiotic_use,inc.parameters_prev,inc.max_year)))
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
        tmp_year = max(ag.cal_year-current_age+asthma_age,min_year)
        while find_asthma_age && asthma_age < 110
            if asthma_age == 3
                has_asthma = rand(Bernoulli(prevalence_equation(tmp_sex,asthma_age,tmp_year,tmp_family_hist,tmp_abx_num,inc.parameters_prev,max_year)))
            else 
                has_asthma = rand(Bernoulli(incidence_equation(tmp_sex,asthma_age,tmp_year,tmp_family_hist,tmp_abx_num,inc.parameters,max_year)))
            end
            if has_asthma
                return asthma_age
            end
            asthma_age += 1 
            asthma_age = min(asthma_age,inc.max_age)
            tmp_year += 1
        end
        return asthma_age 
    end
end



function crude_incidence(sex,age,cal_year,parameters)
    poly_age = poly_age_calculator(age)
    return exp(parameters[:β0] + parameters[:βsex]*sex +  parameters[:βyear]*cal_year + parameters[:βsexyear] * sex * cal_year + sum(parameters[:βage] .* poly_age) + sum(parameters[:βsexage] .* sex .* poly_age))
end

function incidence_equation(sex,age,cal_year,fam_hist,dose,param,max_year)
    correction_year = min(cal_year,max_year+1)
    cal_year = min(cal_year,max_year)
    p0 = crude_incidence(sex,age,cal_year,param)
    return inverse_logit(logit(p0) + fam_hist * log_OR_family_history(age,param[:βfam_hist]) + log_OR_abx_exposure(age,dose,param[:βabx_exp]) + param[:βcorrection][(correction_year,sex,min(age,63))].correction[1])
end

function log_OR_family_history(age,param)
    param[1] + (min(5,age)-3)*param[2]
end

function log_OR_abx_exposure(age,dose,param)
    if (age > 7) | (dose == 0)
        return 0
    else
        return param[1] + param[2] * min(age,7) + param[3] * min(dose,3)
    end
end

function prevalence_equation(sex,age,cal_year,family_hist,dose,param,max_year)
    correction_year = min(cal_year,max_year+1)
    cal_year = min(cal_year,max_year)
    p0 = crude_prevalence(sex,age,cal_year,param)
    return inverse_logit(logit(p0) + family_hist * log_OR_family_history(age,param[:βfam_hist]) + log_OR_abx_exposure(age,dose,param[:βabx_exp]) + param[:βcorrection][(correction_year,sex,min(age,63))].correction[1])
end

function crude_prevalence(sex,age,cal_year,parameters)
    poly_year = poly_year_calculator(cal_year)
    poly_age = poly_age_calculator(age)
    poly_yearage = vec(poly_year .* poly_age')
    return exp(parameters[:β0] + parameters[:βsex]*sex +  sum(parameters[:βyear] .* poly_year) + sum(parameters[:βage] .* poly_age) + sum(parameters[:βsexyear] .* sex .* poly_year)  + sum(parameters[:βsexage] .* sex .* poly_age) +
    sum(parameters[:βyearage] .* poly_yearage) + sum(parameters[:βsexyearage] .* sex .* poly_yearage))
end

function poly_age_calculator(age,alpha=[32.07692,32.42755,32.76123,32.80415,32.54075],nd=[1,520,179636.923076923,47536813.3328764,11589923664.2537,2683688761696.54,594554071731935])
    fs = zeros(6)
    fs[1] =  1/sqrt(nd[2])
    fs[2] = (age-alpha[1]) / sqrt(nd[3])
    for i in 2:5
        fs[i+1] = ((age-alpha[i]) * sqrt(nd[i+1]) * fs[i] - nd[i+1] / sqrt(nd[i]) * fs[i-1]) / sqrt(nd[i+2])
    end
    popfirst!(fs)
    fs
end

function poly_year_calculator(year,alpha=[2009.5,2009.5],nd=[1,520,17290,456456])
    fs = zeros(3)
    fs[1] =  1/sqrt(nd[2])
    fs[2] = (year-alpha[1]) / sqrt(nd[3])
    fs[3] = ((year-alpha[2]) * sqrt(nd[3]) * fs[2] - nd[3] / sqrt(nd[2]) * fs[1]) / sqrt(nd[4])
    popfirst!(fs)
    fs
end

function random_parameter_initialization!(inc::Incidence)
    inc.parameters[:β0] = rand(Normal(inc.hyperparameters[:β0_μ], inc.hyperparameters[:β0_σ]))
end