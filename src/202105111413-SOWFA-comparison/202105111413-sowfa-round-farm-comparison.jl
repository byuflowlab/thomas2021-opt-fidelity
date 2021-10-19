using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using PyPlot; const plt=PyPlot
using DataFrames
using LsqFit
using Colors, ColorSchemes
using Distributions
using CSV

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

function get_data(;journal=false,case="high-ti",opt=false)

    if opt == true
        lesfile = "../../image-generation/image-data/power/turbine-power-$case-les-opt.txt"
        turbine_powers_by_direction_sowfa = zeros((12,38))
        turbine_powers_by_direction_sowfa[:,:] = transpose(readdlm(lesfile, skipstart=1))

    elseif opt == "both"
        lesfile1 = "../../image-generation/image-data/power/turbine-power-$case-les.txt"
        turbine_powers_by_direction_sowfa1 = zeros((12,38))
        turbine_powers_by_direction_sowfa1[:,:] = transpose(readdlm(lesfile1, skipstart=1))

        lesfile2 = "../../image-generation/image-data/power/turbine-power-$case-les-opt.txt"
        turbine_powers_by_direction_sowfa2 = zeros((12,38))
        turbine_powers_by_direction_sowfa2[:,:] = transpose(readdlm(lesfile2, skipstart=1))

        turbine_powers_by_direction_sowfa = [turbine_powers_by_direction_sowfa1, turbine_powers_by_direction_sowfa2]

    else
        lesfile = "../../image-generation/image-data/power/turbine-power-$case-les.txt"
        turbine_powers_by_direction_sowfa = zeros((12,38))
        turbine_powers_by_direction_sowfa[:,:] = transpose(readdlm(lesfile, skipstart=1))
    end
    
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

function run_flow_farm(farmfunc ;x=nothing, y=nothing, use_local_ti=true, nsamplepoints=1, alpha=0.0, verbose=false, windrose="nantucket", shearfirst=true, case="high-ti", ti=0, ws=0, wd=0, opt=false, p=nothing)
    # load FLOWFarm modelset
    if farmfunc === nothing
        filename = "../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$(case)-alldirections.jl"
    
        include(filename)
        farmfunc = wind_farm_setup
    end
    
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = farmfunc(38)
    
    if (x === nothing) && !opt
        turbine_xp = turbine_x .+ 2000.0
        turbine_yp = turbine_y .+ 2000.0
    elseif opt 
        xy = readdlm("../../image-generation/image-data/layouts/opt/optresultsmilestone.csv",',',skipstart=1)
        turbine_xp = xy[:,1]
        turbine_yp = xy[:,2]
        # println("xy", xy)
        # println("turb x: $turbine_xp")
        # println("turb y: $turbine_yp")
    else
        turbine_xp = x 
        turbine_yp = y
    end

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
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints, alpha=alpha, method="sunflower", pradius=1.0, use_perimeter_points=true)

    # set up wind resource 
    if windrose == "nantucket"
        if p !== nothing 
            windspeeds_local = p[:,2]
            ambient_tis_local = p[:,3]
            winddirections_local = winddirections
            windprobabilities_local = windprobabilities
            measurementheight_local = measurementheight
            wind_resource_local = ff.DiscretizedWindResource(winddirections_local, windspeeds_local, windprobabilities_local, measurementheight_local, air_density, ambient_tis_local, wind_shear_model)
        elseif (ti > 0.0) | (ws > 0.0)
            if ti > 0.0
                ambient_tis_local = zeros(length(winddirections)) .+ ti
            else
                ambient_tis_local = ambient_tis
            end
            if ws > 0.0
                windspeeds_local = zeros(length(winddirections)) .+ ws
            else
                windspeeds_local = windspeeds
            end
            if wd != 0.0
                winddirections_local = [winddirections[wd]]
                windspeeds_local = [windspeeds_local[wd]]
                windprobabilities_local = [windprobabilities[wd]]
                measurementheight_local = [measurementheight[wd]]
                ambient_tis_local = [ambient_tis_local[wd]]
            else
                winddirections_local = winddirections
                windprobabilities_local = windprobabilities
                measurementheight_local = measurementheight
            end
            wind_resource_local = ff.DiscretizedWindResource(winddirections_local, windspeeds_local, windprobabilities_local, measurementheight_local, air_density, ambient_tis_local, wind_shear_model)
        else
            wind_resource_local = wind_resource
        end
    else
        windpeed = [8.0]
        winddirection = [270.0*pi/180]
        windprobability = [1.0]
        windheight = [90.0]
        windtis = [ambient_ti]

        # initialize the wind resource definition
        wind_resource_local = ff.DiscretizedWindResource(winddirection, windpeed, windprobability, windheight, air_density, windtis, wind_shear_model)
    end

    
    # run FLOWFarm with local ti
    # println("x: $(turbine_xp) \ny: $(turbine_yp)")
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

function obj_func_internals(points, windspeed, ti, case, wd)
    # get number of points for fitting
    npoints = length(points)

    # get turbine powers with current parameter values
    turbine_powers_flowfarm = run_flow_farm(use_local_ti=true, nsamplepoints=100, ws=windspeed, ti=ti, case=case, wd=wd)
    
    # return just the powers of interest
    powers_to_return = turbine_powers_flowfarm[points]
    # println(powers_to_return)
    return powers_to_return
end

function obj_func_windspeed(points, parameters; ti=0.05589339140297106, case="low-ti", wd=0)
    # println("obj func windspeed $parameters")
    return obj_func_internals(points, parameters[1], ti, case, wd)
end

function obj_func_ti(points, parameters; windspeed=8.466323515677749, case="low-ti", wd=0)
    return obj_func_internals(points, windspeed, parameters[1], case, wd)
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
    # println(fit)

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
    # println(fit)
    tiout = fit.param[1]

    # return wind speed and TI values from tuning
    return wind_speed_out, tiout

end

function state_powers()

    state_powers_ff = zeros(nstates)
    state_powers_ff[:] = sum(turbine_powers_by_direction_ff, dims=2)[:]

end

function errors(first, second; method="absolute", ratedpower=5E6)
    if method == "absolute"
        return first .- second
    elseif method == "normalizedindividually"
        absolute_errors = errors(first, second, method="absolute")
        normalized_errors = absolute_errors./first
        return normalized_errors
    elseif method == "normbyrated"
        absolute_errors = errors(first, second, method="absolute")
        normalized_errors = absolute_errors./ratedpower
        return normalized_errors
    elseif method == "normbyfirst"
        absolute_errors = errors(first, second, method="absolute")
        normalized_errors = absolute_errors./maximum(first)
        return normalized_errors
    elseif method == "normbyrow"
        absolute_errors = errors(first, second, method="absolute")
        normalized_errors = zeros(size(absolute_errors))
        for i in 1:length(first[:,1])
            normalized_errors[i,:] = absolute_errors[i,:]./maximum(absolute_errors[i,:])
        end
        return normalized_errors
    elseif method == "normbyselfrow"
        firstnormed = zeros(size(first))
        secondnormed = zeros(size(second))
        for i in 1:length(first[:,1])
            firstnormed[i,:] = firstnormed[i,:]/maximum(firstnormed[i,:])
            secondnormed[i,:] = secondnormed[i,:]/maximum(secondnormed[i,:])
        end
    
        normalized_errors = errors(firstnormed, secondnormed, method="absolute")
    
        return normalized_errors
    end
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

function custom_color_map()
    colors = [colorant"#BDB8AD", colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    # @pyimport matplotlib.colors as matcolors
    # cmap = matcolors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],"A")

    return plt.ColorMap("BlueGrayOrange", [colors[3],colors[1],colors[4]])
end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1; case="low-ti", tuning="alldirections", opt=false, guess=false, modelsetopt=true)

    # load wind farm information 
    if modelsetopt
        include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-opt.jl")
    else
        include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$case-$tuning.jl")
    end
    # load data
    turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data(case=case, opt=opt)

    if opt == true
        xy = readdlm("../../image-generation/image-data/layouts/opt/optresultsmilestone.csv", ',', skipstart=1)
        println(size(xy[1]))
        turbine_x_local = xy[:,1]
        turbine_y_local = xy[:,2]
        tail = "-opt"
    else
        tail = ""
        turbine_x_local = turbine_x
        turbine_y_local = turbine_y
    end

    if tuning == "alldirections-ws-and-ti-both" || tuning == "all-ws-and-ti-both"
        paramfile = "../202105181144-38-turb-tune-to-sowfa/tuned-parameters-$case-$tuning.csv"
    elseif guess
        paramfile = "../202105181144-38-turb-tune-to-sowfa/tuned-parameters-$case-guess.csv"
    else
        paramfile = "../202105181144-38-turb-tune-to-sowfa/tuned-parameters-$case-$tuning$tail.csv"
    end
    
    p = readdlm(paramfile, ',', skipstart=1)

    # run FLOWFarm
    turbine_powers_by_direction_ff = run_flow_farm(wind_farm_setup, x=turbine_x_local, y=turbine_y_local, use_local_ti=true, 
        nsamplepoints=nsamplepoints, alpha=0.0, verbose=false, windrose="nantucket", shearfirst=true, case=case, p=p)#ti=0.0456610699321765

    # calculate various errors types 
    # absoluteerror = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="absolute")
    normbymax = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyfirst")
    # normbyrated = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyrated")
    # normindividually = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normalizedindividually")
    turberror = normbymax
    data = convert.(Int64, round.(turberror.*100, digits=0))

    # find waked turbines 
    wake_count = ff.wake_count_iec(turbine_x_local, turbine_y_local, winddirections, rotor_diameter)
    dfwc = DataFrame(wake_count',:auto)

    # save data 
    dfff = DataFrame(turbine_powers_by_direction_ff', :auto)

    CSV.write("turbine-power-ff-$(nsamplepoints)pts-$case-$tuning$tail.txt", dfff, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    dfwc = DataFrame(wake_count',:auto)
    CSV.write("turbine-wakes-$case-$tuning$tail.txt", dfwc, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    nturbines = length(turbine_x_local)
    fig, ax = plt.subplots(figsize=(15, 15))
    ticks = minimum(data):5:maximum(data)
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]
    # cs1 = ColorScheme(range(colorant"#0F2080", colorant"#F5793A", length=10)).colors
    # colormap = [(x.r, x.g, x.b) for x in cs1]
    cmap = custom_color_map()
    d = Dict(:shrink => 0.47, :ticks=>ticks, :aspect=>20, :orientation=>"horizontal", :cmap=>cmap)
    
    rowlabels = convert.(Int64, round.((winddirections.*180.0./pi), digits=0))

    im, cbar = heatmap(data, rowlabels, 1:nturbines, ax=ax,
            cbarlabel="Turbine Power Percent Error", cbar_kw=d)

    directionalpowers_ff = reshape(sum(turbine_powers_by_direction_ff,dims=2), 12)
    directionalpowers_sowfa = reshape(sum(turbine_powers_by_direction_sowfa,dims=2), 12)
    directionalerrors = (directionalpowers_sowfa .- directionalpowers_ff)./directionalpowers_sowfa
   
    aep_ff = 365.0.*24.0.*sum(wind_resource.wind_probabilities.*directionalpowers_ff)
    aep_sowfa = 365.0.*24.0.*sum(wind_resource.wind_probabilities.*directionalpowers_sowfa)
    aep_error = (aep_sowfa - aep_ff)/aep_sowfa
    println("AEP SOWFA (GWh): $(round(aep_sowfa*1E-9, digits=2))")
    println("AEP FLOWFarm (GWh): $(round(aep_ff*1E-9, digits=2))")
    println("AEP Error (%): $(round(aep_error.*100,digits=2))")


    fig, ax = plt.subplots()
    ax.bar(wind_resource.wind_directions.*180.0./pi .+ 5, directionalpowers_sowfa.*1E-6, label="SOWFA", width=10)
    ax.bar(wind_resource.wind_directions.*180.0./pi .- 5, directionalpowers_ff.*1E-6, label="FLOWFarm", width=10)
    ax.set(xticks=wind_resource.wind_directions.*180.0./pi, ylim=[0,70])
    plt.xlabel("Direction (deg.)")
    plt.ylabel("Directional Power (MW)")
    plt.legend(frameon=false)

    fig, ax = plt.subplots()
    ax.bar(wind_resource.wind_directions.*180.0./pi, round.(100.0.*directionalerrors, digits=2), width=15)
    ax.set(xticks=wind_resource.wind_directions.*180.0./pi)
    plt.xlabel("Direction (deg.)")
    plt.ylabel("Directional Error (%)");

end

function heatmap(data, row_labels, col_labels; ax=nothing, cbar_kw=Dict(), cbarlabel="", use_cbar=true, labelpixels=true, fontsize=10, vcolor="w", edgecolor="w")
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the pcolormesh call.
    """

    if ax === nothing
        ax = plt.gca()
    end

    # Plot the heatmap
    im = ax.pcolormesh(data, edgecolor=edgecolor, cmap=cbar_kw[:cmap], vmin=-maximum(abs.(cbar_kw[:ticks])), vmax=maximum(abs.(cbar_kw[:ticks])))

    # Create colorbar
    if use_cbar
        cbar = ax.figure.colorbar(im; ax=ax, cbar_kw...)
        cbar.ax.set_xlabel(cbarlabel, rotation=0, va="bottom", labelpad=16)
    else
        cbar = nothing
    end

    # make pixels square
    ax.set(aspect="equal")

    # set tick labels and locations for x axis
    ax.set(xticklabels=col_labels, xticks=(1:length(col_labels)).-0.5)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=true, bottom=false, labeltop=true, labelbottom=false)

    # set tick labels and locations for y axis
    ax.set(yticklabels=row_labels, yticks=(1:length(row_labels)).-0.5)

    # reverse y axis to make the result more table-like
    ax.invert_yaxis()

    # remove ticks
    ax.tick_params(which="minor", top=false, right=false, bottom=false, left=false)
    ax.tick_params(which="major", top=false, right=false, bottom=false, left=false)

    # removes spines
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.spines["left"].set_visible(false)

    # label the pixels
    if labelpixels
        for i = 1:length(row_labels)
            for j = 1:length(col_labels)
                ax.text(j-0.5,i-0.5,data[i,j],
                        ha="center",va="center",
                        size=fontsize,color=vcolor)
            end
        end
    end

    return im, cbar
end