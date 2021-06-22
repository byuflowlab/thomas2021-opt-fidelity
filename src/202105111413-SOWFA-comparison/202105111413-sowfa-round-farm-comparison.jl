using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using Plots
using DataFrames
using LsqFit

function maxk!(ix, a, k; initialized=false)
    partialsortperm!(ix, a, 1:k, rev=true, initialized=initialized)
    @views collect(zip(ix[1:k], a[ix[1:k]]))
end

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

function get_data(;journal=false)
    # load Niayifar LES data 
    lesfile = "../inputfiles/results-thomas-2019/thomas2019-FinalDirectionalGeneratorPowerOutputBaseline.txt"
    sowfa_les_data = readdlm(lesfile, skipstart=0) 
    sowfa_les_data = sowfa_les_data[:,1:5]

    turbine_powers_by_direction_sowfa = format_sowfa_data(sowfa_les_data, 12, 38)
    # println(size(turbine_powers_by_direction_sowfa))
    # load plantenergy data
    if journal
        # confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline_journal.txt" # need to adjust, in kW
        # confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline_journal_updated202106220921.txt"
        # confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline_journal_updated_100rpts.txt"
        # confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline_journal_updated_disNearwake.txt"
        confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline_journal_ti_init_ff_fix.txt"
    else
        confile = "../inputfiles/results-thomas-2019/bp_turb_power_baseline.txt"
    end
    turbine_powers_by_direction_thomas2019 = zeros((12,38))
    turbine_powers_by_direction_thomas2019[:,:] = transpose(readdlm(confile, skipstart=1))# .*1E3 if using first filename

    return turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019
end

function run_flow_farm(;use_local_ti=true, nsamplepoints=1, alpha=0.0, verbose=false, windrose="nantucket", shearfirst=true, filename="../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")
    # load FLOWFarm modelset
    include(filename)
    nturbines = length(turbine_x)
    turbine_xp = turbine_x[1:nturbines] .+ 2000.0
    turbine_yp = turbine_y[1:nturbines] .+ 2000.0
    if verbose
        println("turb x: $turbine_xp")
        println("turb y: $turbine_yp")
    end
    if !use_local_ti
        localtimodel = ff.LocalTIModelNoLocalTI()

        # initialize model set
        model_set_internal = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
    else
        model_set_internal = model_set
    end
    
    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints, alpha=alpha)

    # set up wind resource 
    if windrose == "nantucket"
        wind_resource_local = wind_resource
    else
        windpeed = [8.0]
        winddirection = [270.0*pi/180]
        windprobability = [1.0]
        windheight = [80.0]
        windtis = [ambient_ti]

        # initialize the wind resource definition
        wind_resource_local = ff.DiscretizedWindResource(winddirection, windpeed, windprobability, windheight, air_density, windtis, wind_shear_model)

    end
    # run FLOWFarm with local ti
    turbine_powers_by_direction_ff = ff.calculate_state_turbine_powers(turbine_xp, turbine_yp, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource_local, power_models, model_set_internal,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z, shearfirst=shearfirst)

    return turbine_powers_by_direction_ff

end

# function run_flow_farm(;use_local_ti=true, nsamplepoints=1, windspeedin=8.0, tiin=0.108)

#     # load FLOWFarm modelset
#     include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")

#     # re-set ambient turbulence intensity 
#     ambient_ti = tiin
#     ambient_tis = zeros(nstates) .+ ambient_ti

#     # set wind speed to provided value(s)
#     windspeeds[:] .= windspeedin
#     # initialize the wind resource definition
#     wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

#     # set up local ti model 
#     # if !use_local_ti
#     #     localtimodel = ff.LocalTIModelNoLocalTI()

#     #     # initialize model set
#     #     model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)
#     # end
    
#     # rotor swept area sample points (normalized by rotor radius)
#     rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

#     # run FLOWFarm with local ti
#     turbine_powers_by_direction_ff = ff.calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
#         hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
#         cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
#         rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

#     return turbine_powers_by_direction_ff

# end

function obj_func_internals(points, windspeed, ti)
    # get number of points for fitting
    npoints = length(points)

    # get turbine powers with current parameter values
    turbine_powers_flowfarm = run_flow_farm(use_local_ti=true, nsamplepoints=100, windspeedin=windspeed, tiin=ti)
    
    # return just the powers of interest
    powers_to_return = zeros(npoints)
    for i in 1:npoints
        powers_to_return[i] = turbine_powers_flowfarm[i,Int(points[i])]
    end
    return powers_to_return
end

function obj_func_windspeed(points, parameters; ti=0.05589339140297106)
    return obj_func_internals(points, parameters[1], ti)
end

function obj_func_ti(points, parameters; windspeed=8.466323515677749)
    return obj_func_internals(points, windspeed, parameters[1])
end

function tune_flow_farm()

    # load les data 
    turbine_powers_by_direction_sowfa, _ = get_data()

    # find out how many state
    nstates = 1 #size(turbine_powers_by_direction_sowfa)[1]

    # select front turbines to tune for wind speed 
    front_turbines = zeros(nstates)
    front_turbines_power = zeros(nstates)
    for i in 1:nstates
        front_turbines[i] = Int(argmax(turbine_powers_by_direction_sowfa[i,:]))
        front_turbines_power[i] = turbine_powers_by_direction_sowfa[i,Int(front_turbines[i])]
    end

    # tune front turbines for wind speed 
    fit = curve_fit(obj_func_windspeed, front_turbines, front_turbines_power, [8.0])
    println(fit)

    # select rear turbines to tune for TI 
    rear_turbines = zeros(nstates)
    rear_turbines_power = zeros(nstates)
    for i in 1:nstates
        rear_turbines[i] = Int(argmin(turbine_powers_by_direction_sowfa[i,:]))
        rear_turbines_power[i] = turbine_powers_by_direction_sowfa[i,Int(rear_turbines[i])]
    end
    wind_speed_out = fit.param[1]

    # tune rear turbines for TI 
    fit = curve_fit(obj_func_ti, rear_turbines, rear_turbines_power, [0.108])
    println(fit)
    tiout = fit.param[1]

    # return wind speed and TI values from tuning
    return wind_speed_out, tiout

end

function state_powers()

    state_powers_ff = zeros(nstates)
    state_powers_ff[:] = sum(turbine_powers_by_direction_ff, dims=2)[:]

end

function errors(first, second; method="absolute")
    return first .- second
end

function errors(first, second; method="normbyfirst")
    absolute_errors = errors(first, second, method="absolute")
    normalized_errors = absolute_errors./maximum(first)
    return normalized_errors
end

function errors(first, second; method="normbyrow")
    absolute_errors = errors(first, second, method="absolute")
    normed_errors = zeros(size(absolute_errors))
    for i in 1:length(first[:,1])
        normed_errors[i,:] = absolute_errors[i,:]/maximum(absolute_errors[i,:])
    end
    return normalized_errors
end

function errors(first, second; method="normbyselfrow")
    firstnormed = zeros(size(first))
    secondnormed = zeros(size(second))
    for i in 1:length(first[:,1])
        firstnormed[i,:] = firstnormed[i,:]/maximum(firstnormed[i,:])
        secondnormed[i,:] = secondnormed[i,:]/maximum(secondnormed[i,:])
    end

    normalized_errors = errors(firstnormed, secondnormed, method="absolute")

    return normalized_errors
end

function circleshape(h, k, r)
    theta = LinRange(0, 2*pi, 500)
    return h .+ r*sin.(theta), k.+ r*cos.(theta)
end

function plotturbines!(p, turbinex, turbiney, rotordiameter; linecolor=:black, fillcolor=:blue, markeralpha=0)
    nturbines = length(turbinex)
    for i in 1:nturbines
        plot!(p, circleshape(turbinex[i],turbiney[i],rotordiameter[i]/2.0), seriestype=[:shape],
            linecolor=linecolor, c=fillcolor, legend=false, aspect_ratio=1, alpha=markeralpha)
    end
end

function plot_comparisons(turbinex, turbiney, rotordiameter, comparisons, names)
    nplots = length(comparisons)
    xplots = round(sqrt(nplots))
    yplots = round(sqrt(nplots))
    
    p = plot()

    plotturbines!(p, turbinex, turbiney. rotordiameter)

    display(p)
end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load data
    turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data()

    # run FLOWFarm
    # turbine_powers_by_direction_ffti = run_flow_farm(nsamplepoints=nsamplepoints)

    # calulate various errors types 
    absolute_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019, method="absolute")
    normbyfirst_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019, method="absolute")
    normbyrow_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019, method="normbyrow")
    normbyselfrow_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019, method="normbyselfrow")

    comparison = [absolute_error, normbyfirst_error, normbyrow_error, normbyselfrow_error]
    names = ["absolute", "normbyfirst", "normbyrow", "normbyselfrow"]
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")
    plot_comparisons(turbinex, turbiney, rotordiameter, comparisons, names)

    
    
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