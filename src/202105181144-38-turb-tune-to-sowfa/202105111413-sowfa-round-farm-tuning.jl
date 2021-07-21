using PyPlot: haskey
using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using DataFrames
using CSV
# import PyPlot; const plt = PyPlot

# include functions from other files 
include("../202105111413-SOWFA-comparison/202105111413-sowfa-round-farm-comparison.jl")

function calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::ff.AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, weighted=true)

    nturbines = length(turbine_x)

    wind_probabilities = wind_resource.wind_probabilities

    nstates = length(wind_probabilities)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),typeof(hub_height[1]),typeof(turbine_yaw[1]),
            typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),typeof(rated_power[1]))
    
    turbine_powers_by_direction = zeros(arr_type,(nstates, nturbines))

    for i = 1:nstates

        rot_x, rot_y = ff.rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[i])

        sorted_turbine_index = sortperm(rot_x)

        turbine_velocities = ff.turbine_velocities_one_direction(rot_x, rot_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set, wind_farm_state_id=i, velocity_only=true)

        wt_power = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                            rated_power, rotor_diameter, turbine_velocities, turbine_yaw, wind_resource.air_density, power_models)
        turbine_powers_by_direction[i,:] = wt_power
    end
    
    return turbine_powers_by_direction
end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load Niayifar LES data 
    lesfile = "../inputfiles/results-niayifar-2016/thomas2019-FinalDirectionalGeneratorPowerOutputBaseline.txt"
    sowfa_les_data = readdlm(lesfile, skipstart=0) 
    sowfa_les_data = sowfa_les_data[:,1:5]
    println(size(sowfa_les_data))

    # load FLOWFarm modelset
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")

    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

    # run FLOWFarm with local ti
    turbine_powers_by_direction_ff = calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    state_powers_ff = sum(turbine_powers_by_direction_ff, 2)

    # format sowfa data for plotting
    turbine_powers_by_direction_sowfa = zeros((nstates, nturbines))
    # state_powers_sowfa = sum()
    for i in 1:nstates
        for j in 1:nturbines
            turbine_powers_by_direction_sowfa[i, j] = sowfa_les_data[(i-1)*nturbines + j, 5]
        end
        turbine_powers_by_direction_sowfa[i, :] /= maximum(turbine_powers_by_direction_sowfa[i, :])
        turbine_powers_by_direction_ff[i, :] /= maximum(turbine_powers_by_direction_ff[i,:])
    end
    
    # println(turbine_powers_by_direction_ff)
    # println("BREAK")
    # println(turbine_powers_by_direction_sowfa)
    difference_turbines = turbine_powers_by_direction_sowfa - turbine_powers_by_direction_ff

    heatmap(difference, xticks=1:2:nturbines, yticks=1:nstates, c=:cividis)

end

function find_upstream_turbines(turbinex, turbiney, winddirection, diameter; inverse=false)

    # find wake count for all turbines in given wind direction 
    wake_count = ff.number_of_wakes_iec(turbinex, turbiney, winddirection, diameter)

    if inverse
        # return waked turbines
        return collect(1:length(turbinex))[wake_count .!= 0]
    else
        # return unwaked turbines 
        return collect(1:length(turbinex))[wake_count .== 0]
    end

end

function obj_func_windspeed(points, parameters, inputs::Dict)

    # load inputs
    nsamplepoints = inputs["nsamplepoints"]
    ti = inputs["ti"]
    case = inputs["case"]
    winddirections = inputs["winddirections"]
    upstream_turbines = inputs["upstream_turbines"]

    # get parameter value (wind speed)
    ws = parameters[1]
    println(ws)

    # initialize upstream turbine powers array
    upstream_turbine_powers = []

    # loop through wind directions
    for i = 1:length(winddirections)

        # run flow farm for given direction
        turbine_powers = run_flow_farm(use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws, ti=ti, case=case, wd=i)

        # add result to array 
        push!(upstream_turbine_powers, turbine_powers[upstream_turbines[i]])
    end
    
    # format results for return to LsqFit 
    usp_flat = collect(Iterators.flatten(upstream_turbine_powers))

    if inputs["method"] == "sum"
        return [sum(usp_flat)]
    else
        return usp_flat
    end
end

function obj_func_ti(points, parameters, inputs::Dict)

    # load inputs
    nsamplepoints = inputs["nsamplepoints"]
    ws = inputs["ws"]
    case = inputs["case"]
    winddirections = inputs["winddirections"]
    downstream_turbines = inputs["downstream_turbines"]

    if length(ws) == 1
        ws = ones(length(winddirections)).*ws
    end
    # get parameter value (wind speed)
    ti = parameters[1]
    println("ti = $ti")

    # initialize upstream turbine powers array
    downstream_turbine_powers = []

    # loop through wind directions
    for i = 1:length(winddirections)

        # run flow farm for given direction
        turbine_powers = run_flow_farm(use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws[i], ti=ti, case=case, wd=i)

        # add result to array 
        push!(downstream_turbine_powers, turbine_powers[downstream_turbines[i]])
    end
    
    # format results for return to LsqFit 
    dsp_flat = collect(Iterators.flatten(downstream_turbine_powers))

    if inputs["method"] == "sum"
        println("FF dst_power_sum: $(sum(dsp_flat))")
        return [sum(dsp_flat)]
    else
        return dsp_flat
    end
end

function tune_wind_speed(winddirections, turbine_x, turbine_y, rotor_diameter, 
                            turbine_powers_by_direction_sowfa, case, method)
    # get the number of wind directions and wind turbines 
    ndirections = length(winddirections)

    if method == "directional"
        # initialize arrays for tuned wind speed and TI values
        opt_speeds = zeros(ndirections)

        # tune inflow wind speed for each direction using front turbines
        for i = 1:ndirections
            # find upstream turbines 
            upstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1])
            # println("dir = $(winddirections[i]), upstream turbines = $(upstream_turbines)")

            # get sowfa data for turbines of interest 
            upstream_turbines_power = turbine_powers_by_direction_sowfa[i, upstream_turbines]
            println(upstream_turbines_power)
            # generate appropriate objective function 
            obj_func_windspeed_local(x1, x2) = obj_func_windspeed(x1, x2, ti=0.05589339140297106, case=case, wd=i)

            # run least squared fit to the sowfa data 
            # println("enter objective")
            fit = curve_fit(obj_func_windspeed_local, upstream_turbines, upstream_turbines_power, [8.0])
            
            # store result
            opt_speeds[i] = fit.param[1]

            println("dir: $(round(winddirections[i]*180.0/pi,digits=0)), opt speed: $(opt_speeds[i])")

        end

    elseif method == "all"

        # initialize upstream turbine array
        upstream_turbines = []

        # initalize turbine powers array 
        upstream_powers = []

        # find upstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_upstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1])
            push!(upstream_turbines, dir_upstream_turbines)

            # get and add powers of freestream turbines in this direction
            dir_upstream_powers = turbine_powers_by_direction_sowfa[i, dir_upstream_turbines]
            push!(upstream_powers, dir_upstream_powers)
        end

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ti"=>0.5,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "upstream_turbines"=>upstream_turbines,
                      "method"=>method)

        # initialize correct call signature for objective function 
        obj_func_local(x1, x2) = obj_func_windspeed(x1, x2, inputs)

        # initialize parameters (wind speed)
        p0 = [8.0]

        # reshape data to 1d 
        ust_flat = collect(Iterators.flatten(upstream_turbines))
        usp_flat = collect(Iterators.flatten(upstream_powers))

        # fit the power 
        fit = curve_fit(obj_func_local, ust_flat, usp_flat, p0)

        opt_speeds = [fit.param[1]]

    elseif method == "sum"
        # initialize upstream turbine array
        upstream_turbines = []

        # initalize turbine powers array 
        upstream_powers = []

        # find upstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_upstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1])
            push!(upstream_turbines, dir_upstream_turbines)

            # get and add powers of freestream turbines in this direction
            dir_upstream_powers = turbine_powers_by_direction_sowfa[i, dir_upstream_turbines]
            push!(upstream_powers, dir_upstream_powers)
        end

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                    "ti"=>0.5,
                    "case"=>case,
                    "winddirections"=>winddirections,
                    "upstream_turbines"=>upstream_turbines,
                    "method"=>method)

        # initialize correct call signature for objective function 
        obj_func_local_sum(x1, x2) = obj_func_windspeed(x1, x2, inputs)

        # initialize parameters (wind speed)
        p0 = [8.0]

        # reshape data to 1d 
        ust_power_sum = sum(collect(Iterators.flatten(upstream_powers)))

        # fit the aep 
        fit = curve_fit(obj_func_local_sum, [1], [ust_power_sum], p0)
        opt_speeds = [fit.param[1]]
    end

    return opt_speeds
end

function tune_ambient_ti(winddirections, opt_speeds, turbine_x, turbine_y, rotor_diameter, 
    turbine_powers_by_direction_sowfa, case, method)
    # get the number of wind directions and wind turbines 
    # get the number of wind directions and wind turbines 
    ndirections = length(winddirections)
    opt_tis = 0.0

    if method == "directional"
        # initialize arrays for tuned wind speed and TI values
        opt_tis  = zeros(ndirections)

        # tune inflow wind speed for each direction using front turbines
        for i = 1:ndirections
            # find upstream turbines 
            downstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1], inverse=true)
            # println("dir = $(winddirections[i]), upstream turbines = $(upstream_turbines)")

            # get sowfa data for turbines of interest 
            downstream_turbines_power = turbine_powers_by_direction_sowfa[i, downstream_turbines]
            println(downstream_turbines_power)

            # set inflow wind speed value
            if length(opt_speeds) == 1
                ws = opt_speeds[1]
            else
                ws = opt_speeds[i]
            end

            # generate appropriate objective function 
            obj_func_ti_local(x1, x2) = obj_func_ti(x1, x2, ws=ws, case=case, wd=i)

            # run least squared fit to the sowfa data 
            # println("enter objective")
            fit = curve_fit(obj_func_windspeed_local, downstream_turbines, downstream_turbines_power, [0.07])
            
            # store result
            opt_tis[i] = fit.param[1]

            println("dir: $(round(winddirections[i]*180.0/pi,digits=0)), opt ti: $(opt_tis[i])")

        end

    elseif method == "all"

        # initialize upstream turbine array
        downstream_turbines = []

        # initalize turbine powers array 
        downstream_powers = []

        # find downstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_downstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1], inverse=true)
            push!(downstream_turbines, dir_downstream_turbines)

            # get and store powers of downstream turbines in this direction
            dir_downstream_powers = turbine_powers_by_direction_sowfa[i, dir_downstream_turbines]
            push!(downstream_powers, dir_downstream_powers)
        end

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ws"=>opt_speeds,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "downstream_turbines"=>downstream_turbines,
                      "method"=>method)

        # initialize correct call signature for objective function 
        obj_func_local_all(x1, x2) = obj_func_ti(x1, x2, inputs)

        # initialize parameters (wind speed)
        p0 = [0.07]

        # reshape data to 1d 
        dst_flat = collect(Iterators.flatten(downstream_turbines))
        dsp_flat = collect(Iterators.flatten(downstream_powers))

        # fit the power 
        fit = curve_fit(obj_func_local_all, dst_flat, dsp_flat, p0)

        opt_tis = [fit.param[1]]

    elseif method == "sum"
        # initialize upstream turbine array
        downstream_turbines = []

        # initalize turbine powers array 
        downstream_powers = []

        # find upstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_downstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1])
            push!(downstream_turbines, dir_downstream_turbines)

            # get and add powers of freestream turbines in this direction
            dir_downstream_powers = turbine_powers_by_direction_sowfa[i, dir_downstream_turbines]
            push!(downstream_powers, dir_downstream_powers)
        end

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                        "ws"=>opt_speeds,
                        "case"=>case,
                        "winddirections"=>winddirections,
                        "downstream_turbines"=>downstream_turbines,
                        "method"=>method)

        # initialize correct call signature for objective function 
        obj_func_local_sum(x1, x2) = obj_func_ti(x1, x2, inputs)

        # initialize parameters (wind speed)
        p0 = [0.04]

        # reshape data to 1d 
        dst_power_sum = sum(collect(Iterators.flatten(downstream_powers)))

        println("sowfa power sum: $(dst_power_sum)")

        # fit the aep 
        fit = curve_fit(obj_func_local_sum, [1], [dst_power_sum], p0)

        opt_tis = [fit.param[1]]
    end

    return opt_tis
end

# function tune_ambient_ti(winddirections, opt_speeds, turbine_powers_by_direction_sowfa, case)

#     # get the number of wind directions and wind turbines 
#     ndirections = length(winddirections)
#     opt_tis = zeros(ndirections)

#     # get number of turbines
#     nturbines = length(turbine_x)
    
#     # tune local TI for each direction using all turbines 
#     for i = 1:ndirections 

#         # select all turbines 
#         downstream_turbines = collect(1:nturbines)

#         # get sowfa data for turbines of interest 
#         downstream_turbines_power = turbine_powers_by_direction_sowfa[i, downstream_turbines]

#         # generate appropriate objective function
#         obj_func_ti_local(x1, x2) = obj_func_ti(x1, x2, windspeed=opt_speeds[i], case=case, wd=i)

#         # run least squared fit to the sowfa data 
#         fit = curve_fit(obj_func_ti_local, downstream_turbines, downstream_turbines_power, [0.04])

#         # store result
#         opt_tis[i] = fit.param[1]
#         println("dir: $(round(winddirections[i]*180.0/pi, digits=0)), opt ti: $(opt_tis[i])")

#     end

#     return opt_tis

# end

function tune_flowfarm_to_sowfa(;case="low-ti", method="directional")

    # load sowfa data
    turbine_powers_by_direction_sowfa, _ = get_data(journal=true, case=case)
    # turbine_powers_by_direction_sowfa .*= 1E-3
    # set path to model set 
    model_set_file = "../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$(case).jl"

    # load flowfarm set up
    include(model_set_file)

    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)

    # tune inflow wind speed(s)
    println("tune speeds")
    opt_speeds = tune_wind_speed(winddirections, turbine_x, turbine_y, rotor_diameter, 
                            turbine_powers_by_direction_sowfa, case, method)
    # opt_speeds = [8.055]
    # tune ambient turbulence intensity
    println("tune ti")
    opt_tis = tune_ambient_ti(winddirections, opt_speeds, turbine_x, turbine_y, rotor_diameter, 
                                turbine_powers_by_direction_sowfa, case, method)

    # save results 
    println("results: $opt_speeds $opt_tis")
    if method != "directional"
        opt_speeds = ones(length(winddirections)).*opt_speeds[1]
        opt_tis = ones(length(winddirections)).*opt_tis[1]
    end
    df = DataFrame(dir=(round.(winddirections*180.0./pi, digits=0)), speed=opt_speeds, ti=opt_tis)
    CSV.write("tuned-parameters-$(case).csv", df)

end