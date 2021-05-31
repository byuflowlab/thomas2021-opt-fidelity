using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using Plots
using DataFrames

function format_sowfa_data(sowfa_les_data, nstates, nturbines)

    # initialize sowfa data containers
    turbine_powers_by_direction_sowfa = zeros((nstates, nturbines))

    # put SOWFA data in correct shape
    for i in 1:nstates
        for j in 1:nturbines
            turbine_powers_by_direction_sowfa[i, j] = sowfa_les_data[(i-1)*nturbines + j, 5]
        end
    end

    return turbine_powers_by_direction_sowfa
    
end

function get_data()
    # load Niayifar LES data 
    lesfile = "../inputfiles/results-thomas-2019/thomas2019-FinalDirectionalGeneratorPowerOutputBaseline.txt"
    sowfa_les_data = readdlm(lesfile, skipstart=0) 
    sowfa_les_data = sowfa_les_data[:,1:5]

    turbine_powers_by_direction_sowfa = format_sowfa_data(sowfa_les_data, 12, 38)

    # load plantenergy data
    confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline.txt"
    turbine_powers_by_direction_thomas2019 = zeros((12,38))
    turbine_powers_by_direction_thomas2019[:,:] = transpose(readdlm(confile, skipstart=1)).*1E3

    return turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019
end

function run_flow_farm(;use_local_ti=true, nsamplepoints=1)
    # load FLOWFarm modelset
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")

    if !use_local_ti
        localtimodel = ff.LocalTIModelNoLocalTI()

        # initialize model set
        model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
    end
    
    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

    # run FLOWFarm with local ti
    turbine_powers_by_direction_ff = ff.calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    return turbine_powers_by_direction_ff

end

function state_powers()

    state_powers_ff = zeros(nstates)
    state_powers_ff[:] = sum(turbine_powers_by_direction_ff, dims=2)[:]

end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load data
    turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data()

    # run FLOWFarm
    turbine_powers_by_direction_ffti = run_flow_farm(nsamplepoints=nsamplepoints)
    turbine_powers_by_direction_ffnoti = run_flow_farm(use_local_ti=false,nsamplepoints=nsamplepoints)

    
    # # calculate differences between SOWFA and FlowFarm in each direction
    # differences_directions = (state_powers_sowfa .- state_powers_ff)./state_powers_sowfa
    
    # # print data
    # df = DataFrame(Dir=wind_resource.wind_directions.*180/pi, SOWFA=state_powers_sowfa.*1E-6, FLOWFarm=state_powers_ff.*1E-6, Error=differences_directions.*100)
    # println(df)

    # # calculate and print AEP data 
    # aep_sowfa = sum(state_powers_sowfa.*wind_resource.wind_probabilities*365*24)
    # aep_ff = sum(state_powers_ff.*wind_resource.wind_probabilities*365*24)
    # aep_error = (aep_sowfa - aep_ff)/aep_sowfa
    # println(aep_sowfa.*1E-9, " ", aep_ff.*1E-9, " ", aep_error*100)
    # difference_turbines = (turbine_powers_by_direction_sowfa .- turbine_powers_by_direction_ff)./turbine_powers_by_direction_sowfa
    # colorgrad = cgrad([:red, :white, :blue])
    # println(minimum(difference_turbines))
    # # heatmap(difference_turbines, xticks=1:2:nturbines, yticks=1:nstates, seriescolor=colorgrad, categorical=false)#, clim=(-0.4,0.4))

    # # heatmap of differences between FLOWFarm and PlantEnergy 
    # println(size(thomas2019_bp_data))
    # println(size(turbine_powers_by_direction_ff))
    # diff_pe_ff_turbs = (thomas2019_bp_data - turbine_powers_by_direction_ff)
    # println(size(diff_pe_ff_turbs))
    # heatmap(diff_pe_ff_turbs, xticks=1:2:nturbines, yticks=1:nstates, seriescolor=colorgrad, categorical=false, clim=(-maximum(abs.(diff_pe_ff_turbs)),maximum(abs.(diff_pe_ff_turbs))))
    # # df = DataFrame(Dir=wind_resource.wind_directions.*180/pi, SOWFA=state_powers_sowfa.*1E-6, PlantEnergy=thomas2019_bp_data, FLOWFarm=state_powers_ff.*1E-6, Error=differences_directions.*100)
    # println(df)
end