function Dict_initializer(parameter_names::Union{Nothing,Vector{Symbol}})
    isnothing(parameter_names) ? nothing : Dict(parameter_names .=> missing)
end

function vec_to_dict(v::AbstractArray, ll::AbstractVector)::AbstractDict
    d = Dict()
        for i in eachindex(ll)
            d[ll[i]] = v[i]
        end
    return d
end

function logit(p)
    log.( p ./ (1 .- p) )
end

function inverse_logit(x)
    exp.(x) ./ (1 .+ exp.(x))
end

function set_up(max_age=111,province="BC",starting_year=2001,time_horizon=19,n=100,population_growth_type="LG")
    if province=="BC" || province=="CA"
        
        agent = Agent(false,0,starting_year,1,true,0,false,0,0,nothing,[0,0],[zeros(4),zeros(4)],0,false,false)

        birth = Birth(nothing,nothing)
        @set! birth.estimate = filter([:year, :province, :projection_scenario] => (x, y, z) -> x >= starting_year && y == province && (z == population_growth_type || z == "past"), master_birth_estimate)
        relative(x) = x/birth.estimate.N[1]
        @set! birth.estimate = transform(birth.estimate,:N => relative)
        @set! birth.initial_population =  filter([:year, :province, :projection_scenario] => (x, y, z) -> x == starting_year && y == province && (z == population_growth_type || z == "past"), master_population_initial_distribution)

        death = Death(Dict_initializer([:β0,:β1,:β2]),nothing)
        @set! death.parameters[:β0] =0;
        @set! death.parameters[:β1] =0;
        @set! death.parameters[:β2] =0;
        @set! death.life_table = groupby(select(unstack(
            select(
                select(
                    filter([:year, :province] => (x, y) -> x >= starting_year && y == province, master_life_table),
                Not(:se)),
            Not(:province)),
        :sex,:prob_death),:F,:M,:year),:year)

        emigration = Emigration(nothing, nothing,nothing)
        @set! emigration.table = groupby(select(
            select(
                filter([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year && y == province && z==population_growth_type, master_emigration_table),
            Not(:province)),
        Not(:proj_scenario)),:year)

        immigration = Immigration(nothing, nothing, nothing,nothing,nothing)
        @set! immigration.table = groupby(select(
            select(
                filter([:year, :province, :proj_scenario] => (x, y,z) -> x > starting_year && y == province && z==population_growth_type, master_immigration_table),
            Not(:province)),
        Not(:proj_scenario)),:year)

        incidence = Incidence(Dict_initializer([:β0_μ,:β0_σ]),
        Dict_initializer([:β0,:βsex,:βyear,:βage,:βsexyear,:βsexage,:βfam_hist,:βabx_exp,:βcorrection]),
        Dict_initializer([:β0,:βsex,:βyear,:βage,:βsexyear,:βsexage,:βyearage,:βsexyearage,:βfam_hist,:βabx_exp,:βcorrection]),
        nothing,nothing,nothing)
        
        @set! incidence.hyperparameters[:β0_μ] = 0;
        @set! incidence.hyperparameters[:β0_σ] = 0.00000001;
        @set! incidence.parameters[:β0] = 34.63398846;
        @set! incidence.parameters[:βsex] = -9.52017810;
        @set! incidence.parameters[:βyear] = -0.01967344;
        @set! incidence.parameters[:βage] = [-6.64423331,7.73720625,-5.63121394,3.90920803,-1.39497027];
        @set! incidence.parameters[:βsexyear] = 0.00461397;
        @set! incidence.parameters[:βsexage] = [-4.45607619,4.70483885,-2.61760564,0.79555703,0.95476291];
        @set! incidence.parameters[:βfam_hist] = [log(1.13),0.3619942]
        @set! incidence.parameters[:βabx_exp] = [1.711+0.115, -0.2920745,0.053];
        @set! incidence.parameters[:βcorrection] = groupby(select(filter([:type] => (x) -> x == "inc" ,master_occurrence_correction),Not([:type])),[:year,:sex,:age]);

        @set! incidence.parameters_prev[:β0] = -2.28093577;
        @set! incidence.parameters_prev[:βsex] = -0.10755806;
        @set! incidence.parameters_prev[:βyear] = [2.83586405,-1.18097542];
        @set! incidence.parameters_prev[:βage] = [1.79932480805632, -2.17989374225804, 3.64152189395539, -2.91796538427475, 1.43423653685647];
        @set! incidence.parameters_prev[:βsexyear] = [1.29279956487906, 0.036861276364171];
        @set! incidence.parameters_prev[:βsexage] = [-7.69209530818354, 2.68306716462003, 0.865308192929771, -0.656000992252807, -0.0270826201453694];
        @set! incidence.parameters_prev[:βyearage] = [50.610032709273, 6.51236955045884, -39.4569160874519, 3.69176099747937, 15.9637932343298, -4.79271775804693, -7.14281869955998, 4.18656498490802, -4.88274672641455, -3.3603262281752];
        @set! incidence.parameters_prev[:βsexyearage] = [-3.19896302105009, 7.24422362459046, -25.7979736592919, 0.253623898303176, 11.3848773603672, -2.57625491419054, 7.61284030050534, 4.17111534541718, -15.2128066205219, 3.70514542334455];
        @set! incidence.parameters_prev[:βfam_hist] = [log(1.13),(log(1.13)+log(2.4))/2-log(1.13)];
        @set! incidence.parameters_prev[:βabx_exp] = [1.711+0.115,-0.225,0.053];
        @set! incidence.parameters_prev[:βcorrection] = groupby(select(filter([:type] => (x) -> x == "prev" ,master_occurrence_correction),Not([:type])),[:year,:sex,:age]);
        
        @set! incidence.min_year = collect(keys(incidence.parameters[:βcorrection])[1])[1]+1
        @set! incidence.max_year = collect(keys(incidence.parameters[:βcorrection])[length(incidence.parameters[:βcorrection])])[1] - 1
        @set! incidence.max_age = 63 
        reassessment = Reassessment(nothing)
        @set! reassessment.table =  groupby(filter([:year, :province] => (x, y) -> x >= starting_year && y == province, master_reassessment),:year)

        control = Control(Dict_initializer([:β0_μ,:β0_σ]), Dict_initializer( [:β0,:βage,:βsex,:βsexage,:βsexage2,:βage2, :βDx2,:βDx3,:θ]))
        @set! control.hyperparameters[:β0_μ] = 0;
        @set! control.hyperparameters[:β0_σ] = 1.678728;
        @set! control.parameters[:βage] = 3.5430381;
        @set! control.parameters[:βage2] =-3.4980710;
        @set! control.parameters[:βsexage] = -0.8161495;
        @set! control.parameters[:βsexage2] = -1.1654264;
        @set! control.parameters[:βsex] =  0.2347807;
        @set! control.parameters[:θ] =  [-0.3950; 2.754];

        exacerbation = Exacerbation(Dict_initializer([:β0_μ,:β0_σ]),
        Dict_initializer([:β0,:βage,:βsex,:βasthmaDx,:βprev_exac1,:βprev_exac2,:βcontrol_C,:βcontrol_PC,:βcontrol_UC,:calibration,:min_year]),
        0)
        @set! exacerbation.initial_rate = 0.347;
        @set! exacerbation.hyperparameters[:β0_μ] = 0;
        @set! exacerbation.hyperparameters[:β0_σ] = 0.0000001;
        @set! exacerbation.parameters[:β0_calibration] = 0.0; # 0.056
        @set! exacerbation.parameters[:βage] = 0;
        @set! exacerbation.parameters[:βsex] = 0;
        @set! exacerbation.parameters[:βasthmaDx] = 0;
        @set! exacerbation.parameters[:βprev_exac1] = 0;
        @set! exacerbation.parameters[:βprev_exac2] = 0;
        @set! exacerbation.parameters[:βcontrol_C] =  log(0.1880058);
        @set! exacerbation.parameters[:βcontrol_PC] =  log(0.3760116);
        @set! exacerbation.parameters[:βcontrol_UC] =  log(0.5640174);
        @set! exacerbation.parameters[:βcontrol] =  0;
        @set! exacerbation.parameters[:calibration] = groupby(select(filter([:province] => (x) -> x == province ,exacerbation_calibration),Not([:province])),[:year,:sex]);
        @set! exacerbation.parameters[:min_year] = collect(keys(exacerbation.parameters[:calibration])[1])[1]+1

        exacerbation_severity = Exacerbation_Severity(Dict_initializer([:p0_μ,:p0_σ]), Dict_initializer([:p,:βprev_hosp_ped,:βprev_hosp_adult]))
        @set! exacerbation_severity.hyperparameters[:p0_μ] = [0.495, 0.195, 0.283, 0.026];
        @set! exacerbation_severity.hyperparameters[:p0_σ] = 100;
        @set! exacerbation_severity.parameters[:p] = ones(4)/4;
        @set! exacerbation_severity.parameters[:βprev_hosp_ped] = 1.79
        @set! exacerbation_severity.parameters[:βprev_hosp_adult] = 2.88

        antibioticExposure = AntibioticExposure(Dict_initializer([:β0_μ,:β0_σ]),Dict_initializer([:θ,:β0,:βage,:βsex,:βcal_year,:β2005,:β2005_cal_year,:fix2000,:βfloor,:midtrends]),nothing)
        @set! antibioticExposure.parameters[:θ] =  727.383;
        @set! antibioticExposure.parameters[:β0] = 110.000442;
        @set! antibioticExposure.parameters[:βage] = 0.0;
        @set! antibioticExposure.parameters[:βsex] = 0.249033;
        @set! antibioticExposure.parameters[:βcal_year] = -0.055100;
        @set! antibioticExposure.parameters[:β2005] =  55.033675;
        @set! antibioticExposure.parameters[:β2005_cal_year] =  -0.027437;
        @set! antibioticExposure.parameters[:fixyear] = nothing;
        @set! antibioticExposure.parameters[:βfloor] = 50/1000;
        @set!  antibioticExposure.parameters[:midtrends] = abx_mid_trends;

        familyHistory = FamilyHistory(nothing,Dict_initializer([:p]))
        @set! familyHistory.parameters[:p] = 0.2927242;

        util = Utility(Dict_initializer([:eq5d,:control,:exac]))
        @set! util.parameters[:eq5d] = eq5d
        # disutil
        @set! util.parameters[:control] = [0.06,0.09,0.10]
        # disutil: duration 1 week for mild and two weeks for the rest
        @set! util.parameters[:exac] = [0.32 * 1,  0.44 * 2 , 0.50 * 2 , 0.56 * 2 ] / 52

        cost = Cost(Dict_initializer([:control,:exac]))
        # 1.66 is the exchange rate btw 2018 USD and 2023 CAD Sept
        @set! cost.parameters[:control] = [2372, 2965, 3127]*1.66;
        @set! cost.parameters[:exac] = [130,594, 2425,9900]*1.66;

        sim = Simulation(max_age,province,starting_year,time_horizon,n,population_growth_type,
        agent,
        birth,
        emigration,
        immigration,
        death,
        incidence,
        reassessment,
        control,
        exacerbation,
        exacerbation_severity,
        antibioticExposure,
        familyHistory,
        util,
        cost,
        nothing,
        (;))

        return sim
    else
        error("Province not supported")
    end
end
