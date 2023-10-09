struct Birth <: Birth_Module
    estimate
    trajectory
    initial_population
end

function process(cal_year::Int,cal_year_index::Int,b::Birth)
    return Agent(rand(Bernoulli(b.estimate.prop_male[cal_year_index])),0,cal_year,cal_year_index,true,0,false,nothing,nothing,nothing,[0,0],[zeros(4),zeros(4)],0,false,false)
end

function process(cal_year::Int,cal_year_index::Int,b::Birth,sex,age)
    return Agent(sex,age,cal_year,cal_year_index,true,0,false,nothing,nothing,nothing,[0,0],[zeros(4),zeros(4)],0,false,false)
end

function process_initial(b::Birth,n::Int)
    tmp_n = round.(Int,b.initial_population.prop*n)
    tmp_index = Vector{Int}[]
    for i in eachindex(tmp_n)
        tmp_index = vcat(tmp_index,fill(i,tmp_n[i]))
    end
    tmp_index
    # wsample(1:(nrow(b.initial_population)),b.initial_population.prop,n,replace=true)
end