using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using DataFrames
using Statistics
# using Plots
# pyplot()
using Colors, ColorSchemes
cs1 = ColorScheme([colorant"#BDB8AD",  colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"])

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

    # get state powers at LES locations 
    # initialize the wind resource definition
    nstates_les = length(directions_les)
    windspeeds_les = ones(nstates_les).*8.0
    windprobabilities_les = ones(nstates_les)
    ambient_tis_les = zeros(nstates_les) .+ ambient_ti
    measurementheight_les = zeros(nstates_les) .+ hub_height[1]
    wind_resource_les = ff.DiscretizedWindResource(directions_les, windspeeds_les, windprobabilities_les, measurementheight_les, air_density, ambient_tis_les, wind_shear_model)

    # get state powers at LES directions
    directional_powers_les = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource_les, power_models, model_set,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z, weighted=false)

    # normalize directional power
    normalized_power_ti = directional_powers./max_power_ff
    normalized_power_ti_les = directional_powers_les./max_power_ff

    # # turn off local TI
    local_ti_model = ff.LocalTIModelNoLocalTI()

    # re-initialize model set
    model_set_no_ti = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, local_ti_model)

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

    # get state powers at LES locations 
    directional_powers_no_ti_les = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource_les, power_models, model_set_no_ti,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z, weighted=false)

    # normalize directional power
    normalized_power_no_ti = directional_powers_no_ti./max_power_ff
    normalized_power_no_ti_les = directional_powers_no_ti_les./max_power_ff

    ms = 5
    colors = [colorant"#BDB8AD",  colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    ls = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    lw = 1
    fs = 13

    directions_model_deg = directions_model.*180.0/pi
    directions_les_deg = directions_les.*180.0/pi
    
    open("horns-rev-direction-$nsamplepoints-sample-points.txt", "w") do io
        writedlm(io, ["# (normalized power) model deg, model niayifar, ours no ti, ours with ti"])
        writedlm(io, [directions_model_deg normalized_power_model_niayifar normalized_power_no_ti normalized_power_ti])
    end

    # p = plot(ylim=[0.1, 1.0], legend=:bottomleft, legendcolumns=2, legendfontsize=fs, labelfontsize=fs, tickfontsize=fs, ylabel="Normalized Power", xlabel="Direction", markerstrokewidth=0.0, grid=false, fg_legend=:transparent, background_color=:transparent, foreground_color=:black)
    # scatter!(p, directions_les_deg, normalized_power_les_niayifar, c=colors[1], label="Niayifar 2016 LES", markersize=ms, markershape=:circle, markerstrokewidth=0.0)
    # plot!(p, directions_model_deg, normalized_power_model_niayifar, c=colors[2], label="Niayifar 2016 Model", linestyle=ls[2], linewidth=lw)
    # plot!(p, directions_model_deg, normalized_power_no_ti, c=colors[4], label="FLOWFarm w/o Local TI", linestyle=ls[4], linewidth=lw)
    # plot!(p, directions_model_deg, normalized_power_ti, c=colors[3], label="FLOWFarm w/Local TI", linestyle=ls[3], linewidth=lw)
    # savefig(p, "horns-rev-direction-$nsamplepoints-sample-points.pdf")
    # display(p)

    # get errors 
    dfe = DataFrame(FFNTI=abs.(normalized_power_les_niayifar-normalized_power_no_ti_les), FFTI=abs.(normalized_power_les_niayifar-normalized_power_ti_les))
    println(describe(dfe))

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

    open("horns-rev-rows-$nsamplepoints-sample-points.txt", "w") do io
        writedlm(io, ["# (normalized power) rows normalized_power_les_niayifar normalized_power_model_niayifar normalized_power_averaged_ff_no_ti normalized_power_averaged_ff_ti"])
        writedlm(io, [rows normalized_power_les_niayifar normalized_power_model_niayifar normalized_power_averaged_ff_no_ti normalized_power_averaged_ff_ti])
    end

    ms = 8
    colors = [colorant"#BDB8AD",  colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]

    # fs = 13
    # p = plot(legendfontsize=fs, labelfontsize=fs, tickfontsize=fs, ylabel="Normalized Power", xticks=1:10, xlabel="Row", markerstrokewidth=0.0, grid=false, xlim=[1,10], ylim=[0,1], fg_legend=:transparent, background_color=:transparent, foreground_color=:black)
    # scatter!(p, rows, normalized_power_les_niayifar, markerstrokewidth=0.0, c=colors[1], label="Niayifar 2016 LES", markershape=:circle, markersize=ms)
    # scatter!(p, rows, normalized_power_model_niayifar, markerstrokewidth=0.0, c=colors[2], label="Niayifar 2016 Model", markershape=:utriangle, markersize=ms)
    # scatter!(p, rows, normalized_power_averaged_ff_no_ti, markerstrokewidth=0.0, c=colors[4], label="FLOWFarm w/o Local TI", markershape=:star5, markersize=(ms*1.25))
    # scatter!(p, rows, normalized_power_averaged_ff_ti, markerstrokewidth=0.0, c=colors[3], label="FLOWFarm w/Local TI", markershape=:dtriangle, markersize=ms)
    # savefig(p, "horns-rev-rows-$nsamplepoints-sample-points.pdf")
    # display(p)


    # get errors 
    dfe = DataFrame(FFNTI=abs.(normalized_power_les_niayifar-normalized_power_averaged_ff_no_ti), FFTI=abs.(normalized_power_les_niayifar-normalized_power_averaged_ff_ti))
    println(describe(dfe))
end