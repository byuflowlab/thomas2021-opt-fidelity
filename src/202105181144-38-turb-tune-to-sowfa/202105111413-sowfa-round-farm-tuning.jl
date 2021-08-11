using FLOWFarm: setup_weibull_distribution
using DataFrames: _broadcast_unalias_helper
using PyPlot: haskey
using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using DataFrames
using CSV
using Distributions

using SNOW 
using Snopt 


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

function update_wind_resource!(wind_resource, ws, ti)
    # set up wind resource 
    wind_resource.ambient_tis .= zeros(length(winddirections)) .+ ti
    wind_resource.wind_speeds .= zeros(length(winddirections)) .+ ws
end

function set_up_farm_for_tuning_opt(p0, case, tuning, opt, nsamplepoints)
    if opt == true 
        otxt = "-opt"
    elseif opt == "both"
        otxt = "-both"
    else
        otxt = ""
    end
    filename = "../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$case-$tuning$otxt.jl"
    include(filename)

    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)
    
    if opt != true
        turbine_x .+= 2000.0
        turbine_y .+= 2000.0
    else
        xy = readdlm("../../image-generation/image-data/layouts/opt/optresultsmilestone.csv",',',skipstart=1)
        turbine_x .= xy[:,1]
        turbine_y .= xy[:,2]
    end

    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints, alpha=0.0, method="sunflower", pradius=1.0, use_perimeter_points=true)

    # set up wind resource
    update_wind_resource!(wind_resource, p0[1], p0[2])

    return turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
    rotor_sample_points_y, rotor_sample_points_z

end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load Niayifar LES data 
    lesfile = "../inputfiles/results-niayifar-2016/thomas2019-FinalDirectionalGeneratorPowerOutputBaseline.txt"
    sowfa_les_data = readdlm(lesfile, skipstart=0) 
    sowfa_les_data = sowfa_les_data[:,1:5]
    # println(size(sowfa_les_data))

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

function obj_func_windspeed(points, parameters, inputs::Dict; opt=false)

    # load inputs
    nsamplepoints = inputs["nsamplepoints"]
    ti = inputs["ti"]
    case = inputs["case"]
    winddirections = inputs["winddirections"]
    upstream_turbines = inputs["upstream_turbines"]
    ustp = inputs["ustp"]

    # get parameter value (wind speed)
    ws = parameters[1]
    # println(ws)

    # initialize upstream turbine powers array
    upstream_turbine_powers = []

    # loop through wind directions
    for i = 1:length(winddirections)

        # run flow farm for given direction
        turbine_powers = run_flow_farm(nothing, windrose="nantucket", use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws, ti=ti, case=case, wd=i, opt=opt)

        # add result to array 
        push!(upstream_turbine_powers, turbine_powers[upstream_turbines[i]])
    end
    
    # format results for return to LsqFit 
    if inputs["method"] == "alldirections"
        usp_flat = [sum(upstream_turbine_powers[i]) for i=1:12]
    else
        usp_flat = collect(Iterators.flatten(upstream_turbine_powers))
    end
    
    # println(maximum(usp_flat)," ",maximum(ustp))
    # println(round.((ustp .- usp_flat)./ustp, digits=3))
    ers = (ustp .- usp_flat)
    println("ws: $ws, error-sum: $(sqrt(sum(ers.^2))), error-max: $(maximum(ers))")
    println("difference values: $(ers)")
    if inputs["method"] == "sum"
        return [sum(usp_flat)]
    else
        return usp_flat
    end
end

function obj_func_ti(points, parameters, inputs::Dict; opt=false)

    # load inputs
    nsamplepoints = inputs["nsamplepoints"]
    ws = inputs["ws"]
    case = inputs["case"]
    winddirections = inputs["winddirections"]
    downstream_turbines = inputs["downstream_turbines"]
    dstp = inputs["dstp"]

    if length(ws) == 1
        ws = ones(length(winddirections)).*ws
    end
    # get parameter value (wind speed)
    ti = parameters[1]

    # initialize upstream turbine powers array
    downstream_turbine_powers = []

    # loop through wind directions
    for i = 1:length(winddirections)

        # run flow farm for given direction
        turbine_powers = run_flow_farm(nothing, windrose="nantucket", use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws[i], ti=ti, case=case, wd=i, opt=opt)

        # add result to array 
        push!(downstream_turbine_powers, turbine_powers[downstream_turbines[i]])
    end
    # println("dst: ", downstream_turbine_powers)
    
    # format results for return to LsqFit 
    if inputs["method"] == "alldirections"
        dsp_flat = [sum(downstream_turbine_powers[i]) for i=1:12]
    else
        dsp_flat = collect(Iterators.flatten(downstream_turbine_powers))
    end
    ers = (dstp .- dsp_flat)
    println("ti: $ti, error-sum: $(sqrt(sum(ers.^2))), error-max: $(maximum(abs.(ers)))")
    println("difference values: $(ers)")
    if inputs["method"] == "sum"
        return [sum(dsp_flat)]
    else
        return dsp_flat
    end
end

function obj_func_combined(points, parameters, inputs::Dict; opt=false)

    # load inputs
    nsamplepoints = inputs["nsamplepoints"]
    case = inputs["case"]
    winddirections = inputs["winddirections"]
    powers_sowfa = inputs["powers_flat"]
    powers_sowfa_matrix = inputs["powers_sowfa_matrix"]
    method = inputs["method"]

    # get parameter value (wind speed)
    ws = parameters[1]
    
    # get parameter value (turbulence intensity)
    ti = parameters[2]
    if ti <= 0
        ti = 0.0
    end

    # initialize upstream turbine powers array
    turbine_powers = []
    if opt == "both"
        turbine_powers_opt = []
    end

    # loop through wind directions
    for i = 1:length(winddirections)

        if opt == "both"
            # run flow farm for given direction (base case)
            turbine_powers_temp = run_flow_farm(nothing, windrose="nantucket", use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws, ti=ti, case=case, wd=i, opt=false)
            # add result to array 
            push!(turbine_powers, turbine_powers_temp)

            # run flow farm for given direction (opt case)
            turbine_powers_temp_opt = run_flow_farm(nothing, windrose="nantucket", use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws, ti=ti, case=case, wd=i, opt=true)
            # add result to array 
            push!(turbine_powers_opt, turbine_powers_temp_opt)
        else
            # run flow farm for given direction
            turbine_powers_temp = run_flow_farm(nothing, windrose="nantucket", use_local_ti=true, nsamplepoints=nsamplepoints, ws=ws, ti=ti, case=case, wd=i, opt=opt)
            # add result to array 
            push!(turbine_powers, turbine_powers_temp)
        end
        
    end

    if opt == "both"
        tpm_opt = zeros((12,38))
        for i = 1:12; for j=1:38; tpm_opt[i,j] = turbine_powers_opt[i][j]; end; end
    end

    tpm = zeros((12,38))
    for i = 1:12; for j=1:38; tpm[i,j] = turbine_powers[i][j]; end; end
    
    # format results for return to LsqFit 
    if method == "alldirections" || method == "alldirections-ws-and-ti"
        if opt == "both"
            turbine_powers_flat = append!([sum(tpm[i, :]) for i=1:12],[sum(tpm_opt[i, :]) for i=1:12])
        else
            turbine_powers_flat = [sum(tpm[i, :]) for i=1:12]
        end
        # println(turbine_powers_flat)
        # quit()
        # println(turbine_powers)
        ers = (powers_sowfa .- turbine_powers_flat)
        println("ws: $ws, ti: $ti, sqrt(sum(ers.^2)): $(sqrt(sum(ers.^2))), error-max: $(maximum(ers))")
        # println("difference values: $(ers)")
    else
        # turbine_powers_flat = collect(Iterators.flatten(turbine_powers))
        if opt == "both"
            tpm_flat = append!(collect(Iterators.flatten(tpm)), collect(Iterators.flatten(tpm_opt)))
        else
            tpm_flat = collect(Iterators.flatten(tpm))
        end
        turbine_powers_flat = tpm_flat
        # println("diff: ", (turbine_powers_flat.-tpm_flat))

        # println("ff matrix: ", turbine_powers)
        # println("ff flat: ", size(turbine_powers_flat))
        ers = (powers_sowfa .- tpm_flat)
        # ersmat = powers_sowfa_matrix .- tpm
        println("ws: $ws, ti: $ti, sqrt(sum(ers.^2)): $(sqrt(sum(ers.^2))), error-max: $(maximum(ers))")
        # println("ws: $ws, ti: $ti, sqrt(sum(ersmat.^2)): $(sqrt(sum(ersmat.^2))), error-max: $(maximum(ersmat))")

    end
    
    # println(maximum(usp_flat)," ",maximum(ustp))
    # println(round.((ustp .- usp_flat)./ustp, digits=3))
    
    if inputs["method"] == "sum"
        return [sum(turbine_powers_flat)]
    else
        return turbine_powers_flat
    end
end

function tune_wind_speed(winddirections, turbine_x, turbine_y, rotor_diameter, 
                            turbine_powers_by_direction_sowfa, case, method; opt=false)
    # get the number of wind directions and wind turbines 
    ndirections = length(winddirections)

    if method == "directional"
        # initialize arrays for tuned wind speed and TI values
        opt_speeds = zeros(ndirections)

        # tune inflow wind speed for each direction using front turbines
        for i = 1:ndirections
            # find upstream turbines
            upstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter)
            # println("dir = $(winddirections[i]), upstream turbines = $(upstream_turbines)")

            # get sowfa data for turbines of interest 
            upstream_turbines_power = turbine_powers_by_direction_sowfa[i, upstream_turbines]
            
            # generate appropriate objective function 
            obj_func_windspeed_local(x1, x2) = obj_func_windspeed(x1, x2, ti=0.05589339140297106, case=case, wd=i, opt=opt)

            # run least squared fit to the sowfa data 
            # println("enter objective")
            fit = curve_fit(obj_func_windspeed_local, upstream_turbines, upstream_turbines_power, [8.0])
            
            # store result
            opt_speeds[i] = fit.param[1]

            println("dir: $(round(winddirections[i]*180.0/pi,digits=0)), opt speed: $(opt_speeds[i])")

        end

    elseif method == "all"

        # initialize upstream turbine array
        # upstream_turbines = []

        # initalize turbine powers array 
        upstream_powers = []

        # # find upstream turbines and powers in all wind directions
        upstream_turbines  = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections, rotor_diameter)
        for i = 1:ndirections
            # get and add powers of freestream turbines in this direction
            dir_upstream_powers = turbine_powers_by_direction_sowfa[i, upstream_turbines[i]]
            push!(upstream_powers, dir_upstream_powers)
        end

        # reshape data to 1d 
        ust_flat = collect(Iterators.flatten(upstream_turbines))
        usp_flat = collect(Iterators.flatten(upstream_powers))

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ti"=>0.06,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "upstream_turbines"=>upstream_turbines,
                      "method"=>method,
                      "ustp"=>usp_flat)

        # initialize correct call signature for objective function 
        obj_func_local1(x1, x2) = obj_func_windspeed(x1, x2, inputs, opt=opt)

        # initialize parameters (wind speed)
        p0 = [8.0]

        # fit the power 
        fit = curve_fit(obj_func_local1, ust_flat, usp_flat, p0)

        opt_speeds = [fit.param[1]]

    elseif method == "alldirections"

        # initialize upstream turbine array
        # upstream_turbines = []

        # initalize turbine powers array 
        upstream_powers = []

        # # find upstream turbines and powers in all wind directions
        upstream_turbines  = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections, rotor_diameter)
        for i = 1:ndirections
            # get and add powers of freestream turbines in this direction
            dir_upstream_powers = turbine_powers_by_direction_sowfa[i, upstream_turbines[i]]
            push!(upstream_powers, dir_upstream_powers)
        end

        # reshape data to 1d 
        ust_flat = [sum(upstream_turbines[i]) for i=1:12] 
        usp_flat = [sum(upstream_powers[i]) for i=1:12]
        
        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ti"=>0.06,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "upstream_turbines"=>upstream_turbines,
                      "method"=>method,
                      "ustp"=>usp_flat)

        # initialize correct call signature for objective function 
        obj_func_local2(x1, x2) = obj_func_windspeed(x1, x2, inputs, opt=opt)

        # initialize parameters (wind speed)
        p0 = [8.0]

        # fit the power 
        fit = curve_fit(obj_func_local2, ust_flat, usp_flat, p0)

        opt_speeds = [fit.param[1]]

    elseif method == "sum"
        # initialize upstream turbine array
        upstream_turbines = []

        # initalize turbine powers array 
        upstream_powers = []

        # find upstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_upstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter)
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
        obj_func_local_sum(x1, x2) = obj_func_windspeed(x1, x2, inputs, opt=opt)

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
    turbine_powers_by_direction_sowfa, case, method; opt=false)
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
            downstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1], inverse=true)
            # println("dir = $(winddirections[i]), upstream turbines = $(upstream_turbines)")

            # get sowfa data for turbines of interest 
            downstream_turbines_power = turbine_powers_by_direction_sowfa[i, downstream_turbines]
            # println(downstream_turbines_power)

            # set inflow wind speed value
            if length(opt_speeds) == 1
                ws = opt_speeds[1]
            else
                ws = opt_speeds[i]
            end

            # generate appropriate objective function 
            obj_func_ti_local(x1, x2) = obj_func_ti(x1, x2, ws=ws, case=case, wd=i, opt=opt)

            # run least squared fit to the sowfa data 
            # println("enter objective")
            fit = curve_fit(obj_func_windspeed_local, downstream_turbines, downstream_turbines_power, [0.07])
            
            # store result
            opt_tis[i] = fit.param[1]

            println("dir: $(round(winddirections[i]*180.0/pi,digits=0)), opt ti: $(opt_tis[i])")

        end

    elseif method == "all"

        # initalize turbine powers array 
        downstream_powers = []

        # find downstream turbines and powers in all wind directions
        downstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections, rotor_diameter, inverse=true)
        for i = 1:ndirections
            # get and store powers of downstream turbines in this direction
            dir_downstream_powers = turbine_powers_by_direction_sowfa[i, downstream_turbines[i]]
            push!(downstream_powers, dir_downstream_powers)
        end

        # reshape data to 1d 
        dst_flat = collect(Iterators.flatten(downstream_turbines))
        dsp_flat = collect(Iterators.flatten(downstream_powers))

        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ws"=>opt_speeds,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "downstream_turbines"=>downstream_turbines,
                      "method"=>method,
                      "dstp"=>dsp_flat)

        # initialize correct call signature for objective function 
        obj_func_local_all1(x1, x2) = obj_func_ti(x1, x2, inputs, opt=opt)

        # initialize parameters (wind speed)
        p0 = [0.06]

        # fit the power 
        fit = curve_fit(obj_func_local_all1, dst_flat, dsp_flat, p0)

        opt_tis = [fit.param[1]]

    elseif method == "alldirections"

        # initalize turbine powers array 
        downstream_powers = []

        # find downstream turbines and powers in all wind directions
        downstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections, rotor_diameter, inverse=true)
        for i = 1:ndirections
            # get and store powers of downstream turbines in this direction
            dir_downstream_powers = turbine_powers_by_direction_sowfa[i, downstream_turbines[i]]
            push!(downstream_powers, dir_downstream_powers)
        end

        # reshape data to 1d 
        dst_flat = [sum(downstream_turbines[i]) for i=1:12] #collect(Iterators.flatten(downstream_turbines))
        dsp_flat = [sum(downstream_powers[i]) for i=1:12] #collect(Iterators.flatten(downstream_powers))
        # println("dspsowfa: ", dsp_flat)
        # initialize inputs
        inputs = Dict("nsamplepoints"=>100,
                      "ws"=>opt_speeds,
                      "case"=>case,
                      "winddirections"=>winddirections,
                      "downstream_turbines"=>downstream_turbines,
                      "method"=>method,
                      "dstp"=>dsp_flat)

        # initialize correct call signature for objective function 
        obj_func_local_all2(x1, x2) = obj_func_ti(x1, x2, inputs, opt=opt)

        # initialize parameters (wind speed)
        p0 = [0.06]

        # fit the power 
        fit = curve_fit(obj_func_local_all2, dst_flat, dsp_flat, p0)

        opt_tis = [fit.param[1]]

    elseif method == "sum"
        # initialize upstream turbine array
        downstream_turbines = []

        # initalize turbine powers array 
        downstream_powers = []

        # find upstream turbines and powers in all wind directions
        for i = 1:ndirections
            # get and add freestream turbines in this direction
            dir_downstream_turbines = ff.find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter)
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
        obj_func_local_sum(x1, x2) = obj_func_ti(x1, x2, inputs, opt=opt)

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

function tune_speed_and_ti(winddirections, turbine_x, turbine_y, rotor_diameter, 
                            turbine_powers_by_direction_sowfa, case, method; opt=false,
                            pf=nothing, wt=nothing, nsamplepoints=100)

    # reshape data to 1d 
    if method == "alldirections-ws-and-ti"
        if opt == "both"
            p1 = sum(eachcol(turbine_powers_by_direction_sowfa[1]))
            p2 = sum(eachcol(turbine_powers_by_direction_sowfa[2]))
            powers_flat = append!(p1,p2)
        else
            powers_flat = sum(eachcol(turbine_powers_by_direction_sowfa))
        end
    elseif method == "all-ws-and-ti"
        if opt == "both"
            p1 = collect(Iterators.flatten(turbine_powers_by_direction_sowfa[1]))
            p2 = collect(Iterators.flatten(turbine_powers_by_direction_sowfa[2]))
            powers_flat = append!(p1,p2)
        else
            powers_flat = collect(Iterators.flatten(turbine_powers_by_direction_sowfa))
        end
        
    end
        # println(direction_powers_flat)
    # initialize inputs
    inputs = Dict("nsamplepoints"=>100,
                    "case"=>case,
                    "winddirections"=>winddirections,
                    "powers_flat"=>powers_flat,
                    "method"=>method,
                    "powers_sowfa_matrix"=>turbine_powers_by_direction_sowfa,
                    "pf"=>pf,
                    "wt"=>wt)

    # println("sowfa matrix: ", turbine_powers_by_direction_sowfa)
    println("sowfa flat: ", size(powers_flat))
    # initialize correct call signature for objective function 
    obj_func_local3(x1, x2) = obj_func_combined(x1, x2, inputs, opt=opt)

    # initialize parameters (wind speed)
    p0 = [8.0, 0.07]

    # fit the power 
    if pf === nothing
        fit = curve_fit(obj_func_local3, 1:length(powers_flat), powers_flat, p0)
        opt_speeds = [fit.param[1]]
        opt_tis = [fit.param[2]]
    else # use SNOW and SNOPT to fit 

        turbine_xp, turbine_yp, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y, rotor_sample_points_z = set_up_farm_for_tuning_opt(p0, case, method, opt, nsamplepoints)
        
        inputs["args"] = (turbine_xp, turbine_yp, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y, rotor_sample_points_z)

        inputs["scale"] = 1.0

        function obj_func_snow(g, p, inputs) 
            powers_sowfa = inputs["powers_flat"]
            pf = inputs["pf"]
            wt = inputs["wt"]
            scale = inputs["scale"]
            if p[1] < 0; p[1] = 0; end
            if p[2] < 0; p[2] = 0; end
            (turbine_xp, turbine_yp, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
            rotor_sample_points_y, rotor_sample_points_z) = inputs["args"]

            update_wind_resource!(wind_resource, p[1], p[2])

            turbine_powers_by_direction_ff = ff.calculate_state_turbine_powers(turbine_xp, turbine_yp, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
            rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

            difference = powers_sowfa .- collect(Iterators.flatten(turbine_powers_by_direction_ff))

            # LSE
            f0 = (1/length(powers_sowfa))*sum(abs2, difference)

            # r = sqrt(sum(abs2, wt.*[(p[i] - mean(pf[i]))/std(pf[i]) for i=1:length(p)]))
            r = sum(wt.*[pf[i](p[i]) for i=1:length(p)])

            # println([pf[i](p[i]) for i=1:length(p)])
            # println([(p[i] - mean(pf[i]))/std(pf[i]) for i=1:length(p)])
            # # println(pf)
            # println(p, pf)
            # println((p[1]-mean(pf[1]))/std(pf[1]))
            obj = f0/scale - r

            print(".")
            print(" f0$(f0/scale) r$r obj$obj\n")
            return f0/scale - r
        end
        obj_func_snow(g, p) = obj_func_snow(g, p, inputs) 

        scale = round(obj_func_snow(nothing, p0), digits=0)
        println("scale = $scale")
        inputs["scale"] = scale
        snopt_opt = Dict(
                "Derivative option" => 1,
                "Verify level" => 0,
                "Major optimality tolerance" => 1E-3,
            )

        # initialize solver
        solver = SNOPT(options=snopt_opt)

        # initialize SNOW options
        options = Options(;solver, derivatives=ForwardFD())

        ng = 1
            
        # optimize
        t1 = time()
        xopt, fopt, info, out = minimize(obj_func_snow, p0, ng, -Inf, Inf, -Inf, 0.0, options)
        t2 = time()
        opt_speeds = [xopt[1]]
        opt_tis = [xopt[2]]
        print("\n")
    end

    

    return opt_speeds, opt_tis
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

function tune_flowfarm_to_sowfa(;case="low-ti", method="alldirections-ws-and-ti", opt=false, regular=false, npoints=100)

    # load sowfa data
    turbine_powers_by_direction_sowfa, _ = get_data(journal=true, case=case, opt=opt)
    # turbine_powers_by_direction_sowfa .*= 1E-3
    # set path to model set 
    model_set_file = "../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$case-$method.jl"

    # load flowfarm set up
    include(model_set_file)
    # println(turbine_powers_by_direction_sowfa)
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)

    # if running optimized layout, load and reassign turbine locations
    if opt == true
        xy = readdlm("../../image-generation/image-data/layouts/opt/optresultsmilestone.csv", ',', skipstart=1)
        println(size(xy[1]))
        turbine_x = xy[:,1]
        turbine_y = xy[:,2]
    elseif opt == "both"
        xy = readdlm("../../image-generation/image-data/layouts/opt/optresultsmilestone.csv", ',', skipstart=1)
        turbine_x = [xy[:,1], turbine_x]
        turbine_y = [xy[:,2], turbine_y]
    end

    if regular
        # dws = Normal(8.0, 1.0) # wind speed distribution 
        # dti = Normal(0.06, 0.01) # turbulence intensity distribution
        dws(x) = pdf(Normal(8.0, 1.0), x) # wind speed distribution 
        dti(x) = pdf(Normal(0.06, 0.01), x) # turbulence intensity distribution
        pf = [dws, dti] # array of pdfs for parameters 

        wt = [1.0, 1.0].*1 # array of weights for pdfs
    else
        pf = nothing 
        wt = nothing
    end

    if method == "alldirections-ws-and-ti"
        println("tune speeds and tis")
        opt_speeds, opt_tis = tune_speed_and_ti(winddirections, turbine_x, turbine_y, rotor_diameter, 
        turbine_powers_by_direction_sowfa, case, method, opt=opt, wt=wt, pf=pf, nsamplepoints=100)
    elseif method == "all-ws-and-ti"
        println("tune speeds and tis")
        opt_speeds, opt_tis = tune_speed_and_ti(winddirections, turbine_x, turbine_y, rotor_diameter, 
        turbine_powers_by_direction_sowfa, case, method, opt=opt, wt=wt, pf=pf, nsamplepoints=100)
    else
        # tune inflow wind speed(s)
        println("tune speeds")
        opt_speeds = tune_wind_speed(winddirections, turbine_x, turbine_y, rotor_diameter, 
                                turbine_powers_by_direction_sowfa, case, method, opt=opt)
        # opt_speeds = zeros(12).+8.055
        # opt_speeds = [8.055]
        # tune ambient turbulence intensity
        println("tune ti")
        opt_tis = tune_ambient_ti(winddirections, opt_speeds, turbine_x, turbine_y, rotor_diameter, 
                                    turbine_powers_by_direction_sowfa, case, method, opt=opt)
    end
    # save results 
    println("results: $opt_speeds $opt_tis")
    if method != "directional"
        opt_speeds = ones(length(winddirections)).*opt_speeds[1]
        opt_tis = ones(length(winddirections)).*opt_tis[1]
    end
    df = DataFrame(dir=(round.(winddirections*180.0./pi, digits=0)), speed=opt_speeds, ti=opt_tis)
    if regular
        reg = "-reg" 
    else
        reg = "" 
    end
    if opt == true
        CSV.write("tuned-parameters-$(case)-$(method)-opt$reg.csv", df)
    elseif opt == "both"
        CSV.write("tuned-parameters-$(case)-$(method)-$opt$reg.csv", df)
    else
        CSV.write("tuned-parameters-$(case)-$(method)$reg.csv", df)
    end
end