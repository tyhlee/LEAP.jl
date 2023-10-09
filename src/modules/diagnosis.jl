struct Diagnosis <: Diagnosis_Module
    table
    table_mis
end

function process(ag::Agent,diag::Diagnosis)
    max_year = length(diag.table)
    if ag.age < 4
        ag.has_asthma
    else 
        if ag.has_asthma # if yes asthma, can get mis dx of not having asthma
            rand(Bernoulli(diag.table[min(max_year,ag.cal_year_index)][ag.age-3,ag.sex+3]))
        else # if no asthma but then can get mis dx of asthma
            rand(Bernoulli(diag.table_mis[min(ag.cal_year_index,max_year)][ag.age-3,ag.sex+3]))
        end
    end
end