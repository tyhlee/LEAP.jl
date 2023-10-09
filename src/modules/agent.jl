# better to use UInt64 but the syntax may not be intuitive to others
struct Agent  <: Agent_Module
    sex::Bool # true: male
    age::Int
    cal_year::Int
    cal_year_index::Int
    alive::Bool # true: alive
    num_antibiotic_use::Int
    has_asthma::Bool # true: has asthma
    asthma_age::Union{Nothing,Int} # asthma dx age
    severity::Union{Nothing,Int} # asthma severity level: 1 (mild), 2 (severe), 3 very severe
    control::Union{Nothing,Vector{Float64}} # asthma control level: 1 (controlled), 2 (partially controlled), 3 (fully controlled)
    exac_hist::Union{Nothing,Vector{Int}} # total number of exacerbations
    exac_sev_hist::Union{Nothing,Vector{Vector{Int}}} # number of exacerbations by severity
    total_hosp::Int # total number of very severe asthma exacerbation "hospitalization" 
    family_hist::Bool
    asthma_status::Bool
end

function set_agent!(ag, sex, age,cal_year,cal_year_index,alive, num_antibiotic_use,
    has_asthma,asthma_age,asthma_severity,asthma_control,asthma_exac_hist,asthma_exac_sev_history,total_hosp,fam_hist,asthma_status)
    ag.sex = sex
    ag.age= age
    ag.cal_year = cal_year
    ag.cal_year_index = cal_year_index
    ag.alive = alive
    ag.num_antibiotic_use=num_antibiotic_use
    ag.has_asthma = has_asthma
    ag.asthma_age = asthma_age
    ag.severity = asthma_severity
    ag.control = asthma_control
    ag.exac_hist = asthma_exac_hist
    ag.exac_sev_hist = asthma_exac_sev_history
    ag.total_hosp = total_hosp
    ag.family_hist = fam_hist
    ag.asthma_status = asthma_status
    nothing
end

function process_initial(ag::Agent,asthma_age_data)
    if ag.age==0
        return 0
    else
        return StatsBase.sample(Weights(asthma_age_data[1:ag.age+1,Int(ag.sex)+1]))-1
    end
end
