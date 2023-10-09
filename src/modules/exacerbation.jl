struct Exacerbation <: Exacerbation_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
    initial_rate::Float64
end

function process(ag::Agent,exac::Exacerbation)
    parameters = exac.parameters
    tmp_cal_year = max(ag.cal_year,exac.parameters[:min_year])
    tmp_age = min(ag.age,90)
    exacerbation_prediction(parameters[:β0]+
    parameters[:β0_calibration] +
    ag.age * parameters[:βage]+
    ag.sex * parameters[:βsex]+
    ag.control[1] * parameters[:βcontrol_C] +
    ag.control[2] * parameters[:βcontrol_PC] +
    ag.control[3] * parameters[:βcontrol_UC] + log(parameters[:calibration][(tmp_cal_year,Int(ag.sex))][tmp_age-2,"calibrator_multiplier"])) 
end

function process_initial(ag::Agent,exac::Exacerbation)
    parameters = exac.parameters
    tmp_year = max(exac.parameters[:min_year],ag.cal_year-1)
    tmp_age = min(ag.age-1,90)
    if tmp_age <3
        return 0 
    else
        return exacerbation_prediction(parameters[:β0]+
        tmp_age * parameters[:βage]+
        ag.sex * parameters[:βsex]+
        ag.control[1] * parameters[:βcontrol_C] +
        ag.control[2] * parameters[:βcontrol_PC] +
        ag.control[3] * parameters[:βcontrol_UC] + log(parameters[:calibration][(tmp_year,Int(ag.sex))][tmp_age-2,"calibrator_multiplier"]))
    end
end

function random_parameter_initialization!(exac::Exacerbation)
    exac.parameters[:β0] = rand(Normal(exac.hyperparameters[:β0_μ], exac.hyperparameters[:β0_σ]))
end

function exacerbation_prediction(eta::Float64;inv_link::Function=exp)
    rand(Poisson(inv_link(eta)))
end