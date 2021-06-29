using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
import PyPlot; const plt=PyPlot
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

function run_flow_farm(;use_local_ti=true, nsamplepoints=1, alpha=0.0, verbose=false, windrose="nantucket", shearfirst=true, filename="../inputfiles/model-sets/round-farm-38-turbs-12-dirs-high-ti.jl")
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

function custum_color_map()
    colors = [colorant"#BDB8AD", colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    # @pyimport matplotlib.colors as matcolors
    # cmap = matcolors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],"A")

    return plt.ColorMap("BlueGrayOrange", [colors[3],colors[1],colors[4]])
end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load wind farm information 
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-high-ti.jl")

    # load data
    turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data()

    # run FLOWFarm
    turbine_powers_by_direction_ff = run_flow_farm(nsamplepoints=nsamplepoints)

    # calulate various errors types 
    absoluteerror = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="absolute")
    normbymax = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyfirst")
    normbyrated = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyrated")
    normindividually = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normalizedindividually")

    data = convert.(Int64, round.(normbymax.*100, digits=0))

    # save data 
    dfff = DataFrame(turbine_powers_by_direction_ff', :auto)
    CSV.write("turbine_power_ff_$(nsamplepoints)pts.txt", dfff, header=string.(round.(winddirections.*180.0./pi, digits=0)))

    nturbines = length(turbine_x)
    fig, ax = plt.subplots(figsize=(15, 15))
    ticks = minimum(data):5:maximum(data)
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]
    # cs1 = ColorScheme(range(colorant"#0F2080", colorant"#F5793A", length=10)).colors
    # colormap = [(x.r, x.g, x.b) for x in cs1]
    cmap = custum_color_map()
    d = Dict(:shrink => 0.47, :ticks=>ticks, :aspect=>20, :orientation=>"horizontal", :cmap=>cmap)
    
    rowlabels = convert.(Int64, round.((winddirections.*180.0./pi), digits=0))

    im, cbar = heatmap(data, rowlabels, 1:nturbines, ax=ax,
            cbarlabel="Turbine Power Percent Error", cbar_kw=d)
    
    # cbar.set_label("Turbine Power Error, $\%$", rotation=90)
    # ax.set_ylabel("Direction, degrees")
    # ax.set_xlabel("Turbine")
    # plot error on wind farm 
    # fig, ax = plt.subplots()
    # ff.plotlayout!(ax, turbine_x, turbine_y, rotor_diameter)
    # ax.set(xlim=[-2500, 2500], ylim=[-2500, 2500], aspect="equal")

    # fig, ax = plt.subplots(figsize=(15, 15))
    # im, cbar = heatmap(turb_error, 1:nturbines, wind_directions, ax=ax, cmap="bwr", cbarlabel="Turbine Power Error")
    # # cbar.set_label("Turbine Power Error, $\%$", rotation=90)
    # ax.set_ylabel("Direction, degrees")
    # # ax.set_xlabel("Turbine")
    # fig, ax = plt.subplots()
    # im = ax.imshow(norbed_by_max)

    # # We want to show all ticks...
    # ax.set_yticks(1:length(winddirections))
    # ax.set_xticks(1:nturbines)
    # # ... and label them with the respective list entries
    # ax.set_yticklabels(winddirections)
    # ax.set_xticklabels(1:nturbines)

    # # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
    #         rotation_mode="anchor")

    # # Loop over data dimensions and create text annotations.
    # for i = 1:nturbines
    #     for j = 1:length(winddirections)
    #         text = ax.text(j, i, norbed_by_max[i, j],
    #                     ha="center", va="center", color="w")



    #     end
    # end

    # ax.set_title("Harvest of local farmers (in tons/year)")
    # fig.tight_layout()
    # plt.show()

    # normalizedindividually = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normalizedindividually")
    # normbyrow_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyrow")
    # normbyselfrow_error = errors(turbine_powers_by_direction_sowfa, turbine_powers_by_direction_ff, method="normbyselfrow")

    # comparison = [absolute_error, normbyfirst_error, normbyrow_error, normbyselfrow_error]
    # names = ["absolute", "normbyfirst", "normbyrow", "normbyselfrow"]
    # include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")
    # plot_comparisons(turbinex, turbiney, rotordiameter, comparisons, names)

    
    
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
    im = ax.pcolormesh(data, edgecolor=edgecolor, cmap=cbar_kw[:cmap], vmin=minimum(cbar_kw[:ticks]), vmax=maximum(cbar_kw[:ticks]))

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

# function heatmap(data, row_labels, col_labels; ax=nothing, cbar_kw=Dict(), cbarlabel="", use_cbar=true)
#     """
#     Create a heatmap from a numpy array and two lists of labels.

#     Arguments:
#         data       : A 2D numpy array of shape (N,M)
#         row_labels : A list or array of length N with the labels
#                      for the rows
#         col_labels : A list or array of length M with the labels
#                      for the columns
#     Optional arguments:
#         ax         : A matplotlib.axes.Axes instance to which the heatmap
#                      is plotted. If not provided, use current axes or
#                      create a new one.
#         cbar_kw    : A dictionary with arguments to
#                      :meth:`matplotlib.Figure.colorbar`.
#         cbarlabel  : The label for the colorbar
#     All other arguments are directly passed on to the imshow call.
#     """
#     df = DataFrame(data, :auto)
#     println(df)
#     if ax === nothing
#         ax = plt.gca()
#     end

#     # Plot the heatmap
#     im = ax.imshow(data, interpolation="none")

#     # Create colorbar
#     if use_cbar
#         cbar = ax.figure.colorbar(im; ax=ax, cbar_kw...)
#         cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
#     else
#         cbar = nothing
#     end
#     # We want to show all ticks...
#     ax.set_xticks(0:(size(data)[2]-1))
#     ax.set_yticks(0:(size(data)[1]-1))
#     # ... and label them with the respective list entries.
#     println(col_labels)
#     println(row_labels)
#     ax.set_xticklabels(col_labels)
#     ax.set_yticklabels(row_labels)

#     # Let the horizontal axes labeling appear on top.
#     ax.tick_params(top=true, bottom=false, labeltop=true, labelbottom=false)

#     # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=-0, ha="center", rotation_mode="anchor")

#     # Turn spines off and create white grid.
#     ax.spines["right"].set_visible(false)
#     ax.spines["top"].set_visible(false)
#     ax.spines["bottom"].set_visible(false)
#     ax.spines["left"].set_visible(false)

#     ax.set_xticks((1:size(data)[2]).-.5, minor=true)
#     ax.set_yticks((1:size(data)[1]).-.5, minor=true)
#     ax.grid(which="minor", color="w", linestyle="-", linewidth=3)

#     ax.tick_params(which="minor", top=false, right=false, bottom=false, left=false)
#     ax.tick_params(which="major", top=false, right=false, bottom=false, left=false)

#     return im, cbar
# end