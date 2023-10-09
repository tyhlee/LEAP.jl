struct Control <: Control_Module
    hyperparameters::Union{AbstractDict,Nothing}
    parameters::Union{AbstractDict,Nothing}
end

function process(ag::Agent,ctl::Control)
    age_scaled = ag.age / 100
    control_prediction(ctl.parameters[:β0]+
    age_scaled*ctl.parameters[:βage]+
    ag.sex*ctl.parameters[:βsex]+
    age_scaled * ag.sex * ctl.parameters[:βsexage] +
    age_scaled^2 * ag.sex * ctl.parameters[:βsexage2] +
    age_scaled^2 * ctl.parameters[:βage2], ctl.parameters[:θ])
end


function process_initial(ag::Agent, ctl::Control)
    age_scaled = (ag.age-1) / 100
    control_prediction(ctl.parameters[:β0]+
    age_scaled*ctl.parameters[:βage]+
    ag.sex*ctl.parameters[:βsex]+
    age_scaled * ag.sex * ctl.parameters[:βsexage] +
    age_scaled^2 * ag.sex * ctl.parameters[:βsexage2] +
    age_scaled^2 * ctl.parameters[:βage2], ctl.parameters[:θ])
end

# pred function
function control_prediction(eta::Float64,theta::Union{Float64,Vector{Float64}};inv_link::Function=StatsFuns.logistic)::Union{Float64,Vector{Float64}}
    theta = [-1e5;theta;1e5]
    [inv_link(theta[j+1] - eta) - inv_link(theta[j] - eta) for j in 1:(length(theta)-1)]
end

# input: initialized hyperparameters, empty parameters
function random_parameter_initialization!(ctl::Control)
    ctl.parameters[:β0] = rand(Normal(ctl.hyperparameters[:β0_μ],ctl.hyperparameters[:β0_σ]))
end