struct Reassessment <: Reassessment_Module
    table
end

function process(ag::Agent,ra::Reassessment)
    max_year = length(ra.table)
    if ag.age <4
        ag.has_asthma
    else
        rand(Bernoulli(ra.table[min(ag.cal_year_index,max_year)][ag.age-3,ag.sex+3]))
    end
end