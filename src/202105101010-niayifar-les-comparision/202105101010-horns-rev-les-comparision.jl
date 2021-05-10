using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using Plots

# function to compare directions 
function horns_rev_directions(nsamplepoints=1)

    # load Niayifar LES data 
    lesfile = "../inputfiles/results-niayifar-2016/horns-rev-normalized-power-by-direction-les.txt"
    niayifar_les_data = readdlm(lesfile, ',', skipstart=1)
    normalized_power_les_niayifar = niayifar_les_data[:,2]
    directions_les = niayifar_les_data[:,1].*pi./180.0
    ndirections_les = length(directions_les)
    
    # load Niayifar model data
    modelfile = "../inputfiles/results-niayifar-2016/horns-rev-normalized-power-by-direction-model.txt"
    niayifar_model_data = readdlm(modelfile,  ',', skipstart=1)
    normalized_power_model_niayifar = niayifar_model_data[:, 2]
    directions_model = niayifar_model_data[:,1].*pi./180.0
    ndirections_model = length(directions_model)
    

    # load FLOWFarm modelset
    include("../inputfiles/model-sets/horns-rev-by-direction.jl")

    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

    # run FLOWFarm with local ti
    # get max turbine power
    turbine_inflow_velcities = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set, wind_farm_state_id=1)

    turbinepowers = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, 
        rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, turbine_yaw, air_density, 
        power_models)

    max_power_ff = maximum(turbinepowers).*nturbines

    # get state powers
    directional_powers = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z, weighted=false)

    # normalize directional power
    normalized_power_ti = directional_powers./max_power_ff

    # # turn off local TI
    local_ti_model = ff.LocalTIModelNoLocalTI()

    # re-initialize model set
    model_set_no_ti = ff.WindFarmModelSet(wake_deficit_model, wake_deflection_model, wake_combination_model, local_ti_model)

    # run FLOWFarm without local ti
    # get max turbine power
    turbine_inflow_velcities = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set_no_ti, wind_farm_state_id=1)

    turbinepowers = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, 
        rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, turbine_yaw, air_density, 
        power_models)

    max_power_ff = maximum(turbinepowers).*nturbines

    # get state powers
    directional_powers_no_ti = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set_no_ti,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z, weighted=false)

    # normalize directional power
    normalized_power_no_ti = directional_powers_no_ti./max_power_ff
    scatter(directions_les, normalized_power_les_niayifar, label="LES",legend=:bottomright)
    scatter!(directions_model, normalized_power_model_niayifar, label="Niayifar 2016")
    scatter!(winddirections, normalized_power_ti, label="FLOWFarm w/TI")
    scatter!(winddirections, normalized_power_no_ti, label="FLOWFarm w/o TI", xlabel="Direction", ylabel="Normalized Power",fg_legend=:transparent, title="Rotor Points: $nsamplepoints")
    
end

# function to compare rows
function horns_rev_rows(nsamplepoints=1)

    # load Niayifar LES data 
    lesfile = "../inputfiles/results-niayifar-2016/horns-rev-normalized-power-by-row-les.txt"
    niayifar_les_data = readdlm(lesfile, ',', skipstart=1)
    normalized_power_les_niayifar = niayifar_les_data[:,2]
    rows = niayifar_les_data[:,1]
    nrows = length(rows)
    
    # load Niayifar model data
    modelfile = "../inputfiles/results-niayifar-2016/horns-rev-normalized-power-by-row-model.txt"
    niayifar_model_data = readdlm(modelfile,  ',', skipstart=1)
    normalized_power_model_niayifar = niayifar_model_data[:, 2]

    # load FLOWFarm modelset
    include("../inputfiles/model-sets/horns-rev-by-row.jl")

    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

    # run FLOWFarm with local ti
    turbine_inflow_velcities = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set)

    turbinepowers = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, 
        rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, turbine_yaw, air_density, 
        power_models)
    
    normalizedpower = turbinepowers./maximum(turbinepowers)
    normalized_power_averaged_ff_ti = zeros(nrows)
    for row in 1:nrows
        normalized_power_averaged_ff_ti[row] = mean([normalizedpower[(10 * 0 + 40) + row],
                                                   normalizedpower[(10 * 1 + 40) + row],
                                                   normalizedpower[(10 * 2 + 40) + row]])
    end

    # # turn off local TI
    local_ti_model = ff.LocalTIModelNoLocalTI()

    # re-initialize model set
    model_set_no_ti = ff.WindFarmModelSet(wake_deficit_model, wake_deflection_model, wake_combination_model, local_ti_model)

    # run FLOWFarm with local ti
    
    turbine_inflow_velcities = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set_no_ti)

    turbinepowers = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, 
        rated_speed, rated_power, rotor_diameter, turbine_inflow_velcities, turbine_yaw, air_density, 
        power_models)
    
    normalizedpower = turbinepowers./maximum(turbinepowers)
    normalized_power_averaged_ff_no_ti = zeros(nrows)
    for row in 1:nrows
        normalized_power_averaged_ff_no_ti[row] = mean([normalizedpower[(10 * 0 + 40) + row],
                                                   normalizedpower[(10 * 1 + 40) + row],
                                                   normalizedpower[(10 * 2 + 40) + row]])
    end

    scatter(rows, normalized_power_les_niayifar, label="LES", xlim=[1,10], title="Rotor Points: $nsamplepoints")
    scatter!(rows, normalized_power_model_niayifar, label="Niayifar 2016")
    scatter!(rows, normalized_power_averaged_ff_ti, label="FLOWFarm w/TI")
    scatter!(rows, normalized_power_averaged_ff_no_ti, label="FLOWFarm w/o TI", xlabel="Row", xticks=1:10, ylabel="Normalized Power",fg_legend=:transparent, xlim=[1,10], title="Rotor Points: $nsamplepoints")
    
end