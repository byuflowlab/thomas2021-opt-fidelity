using CSV: length
using SNOW
using Snopt
using DelimitedFiles 
using PyPlot
using DataFrames
using CSV
# using Distributed
using ProgressMeter
using StatsBase

# import model set with wind farm and related details
include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-opt.jl")

# set globals struct for use in wrapper functions
mutable struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_center
    boundary_radius
    obj_scale
    con_scale_boundary
    xyscale
    hub_height
    turbine_yaw
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    wind_resource
    power_models
end

function get_sample_size()
    mean = 13.280
    std = 0.341
    delta = 0.01
    alpha = 0.05 
    beta = 0.05
    z = zscore(2.0)
    println(z)
    n = z*std/del
    return n
end

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    boundary_center = params.boundary_center
    boundary_radius = params.boundary_radius
    con_scale_boundary = params.con_scale_boundary

    # find the number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundary_center, boundary_radius, turbine_x, turbine_y).*con_scale_boundary
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant globals
    rotor_diameter = params.rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
function aep_wrapper(x, params)
    # include relevant globals
    turbine_z = params.turbine_z
    rotor_diameter = params.rotor_diameter
    hub_height = params.hub_height
    turbine_yaw =params.turbine_yaw
    ct_models = params.ct_models
    generator_efficiency = params.generator_efficiency
    cut_in_speed = params.cut_in_speed
    cut_out_speed = params.cut_out_speed
    rated_speed = params.rated_speed
    rated_power = params.rated_power
    wind_resource = params.wind_resource
    power_models = params.power_models
    model_set = params.model_set
    rotor_points_y = params.rotor_points_y
    rotor_points_z = params.rotor_points_z
    obj_scale = params.obj_scale

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as an array
    return AEP
end

# set up optimization problem wrapper function
function wind_farm_opt!(g, x, params; xhistory=nothing)

    nturbines = Int(length(x)/2)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x, params)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x, params)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[1:(end-nturbines)] = spacing_con[:]
    g[end-nturbines+1:end] = boundary_con[:]
    
    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x, params)[1]
    # println("   $AEP")
    # flush(stdout)
    # save steps to history
    if xhistory !== nothing 
        xfloat = [x[i].value for i in 1:length(x)]
        push!(xhistory, xfloat)
    end
    
    return AEP #, dAEP_dx, dcdx, fail
end

function set_up_base_params(params; nrotorpoints=100, alpha=0)

    params_base = deepcopy(params)

    # set sample points 
    rotor_points_y, rotor_points_z = ff.rotor_sample_points(nrotorpoints, method="sunflower", pradius=1, alpha=alpha)

    params_base.rotor_points_y = rotor_points_y
    params_base.rotor_points_z = rotor_points_z

    return params_base

end

function turbine_power_by_direction_from_params(turbine_x, turbine_y, params)
    turbine_powers_by_direction = ff.calculate_state_turbine_powers(turbine_x, turbine_y, params.turbine_z, params.rotor_diameter,
        params.hub_height, params.turbine_yaw, params.ct_models, params.generator_efficiency, params.cut_in_speed,
        params.cut_out_speed, params.rated_speed, params.rated_power, params.wind_resource, params.power_models, params.model_set,
        rotor_sample_points_y=params.rotor_points_y, rotor_sample_points_z=params.rotor_points_z)
    return turbine_powers_by_direction
end

function rerun(layoutid; case="high-ti", tuning="sowfa-nrel", plotresults=false, verbose=true, wec=true, nrotorpoints=1, alpha=0, outdir="./", layoutdir="../inputfiles/farms/startinglayouts/individual/", lspacing=3.0, n="4")

    # get wind farm setup
    diam, turbine_x_start, turbine_y_start, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, _, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38, case=case, tuning=tuning, layoutid=layoutid, nrotorpoints=nrotorpoints, alpha=alpha, layoutdir=layoutdir, lspacing=lspacing)

    # get starting layout
    _, turbine_x_base, turbine_y_base, _, _, _, _, _, 
    _, _, _, _, _, 
    _, _, _, _, _, 
    _, _, _, _, _, _, 
    _, _, _, _, _, 
    _, _, _, _ = wind_farm_setup(38, case=case, tuning=tuning, layoutid=1, nrotorpoints=nrotorpoints, alpha=alpha, layoutdir=layoutdir, lspacing=lspacing)

    optlayoutfile = "../../image-generation/image-data/layouts/opt/$case-$tuning-opt$n-layout$layoutid-aec-wec.csv"

    # extract data
    df = DataFrame(CSV.File(optlayoutfile, skipto=2, header=false))

    # name columns 
    rename!(df,:Column1 => :x,:Column2 => :y)

    turbine_x_opt = df.x
    turbine_y_opt = df.y

    # get number of turbines 
    nturbines = length(turbine_x)

    # scale objective to be between 0 and 1
    obj_scale = 1E-9 #1E-11
    con_scale_boundary = 1E-4
    xyscale = 1 #1E4

    # set wind farm boundary parameters
    boundary_center = [0.0,0.0]
    boundary_radius = 0.5 * (4000. - rotor_diameter[1])  # 1936.8 
    if verbose
        println("boundary radius: $(boundary_radius)")
    end

    params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
        rotor_diameter, boundary_center, boundary_radius, obj_scale, con_scale_boundary, xyscale, hub_height, turbine_yaw, 
        ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
        wind_resource, power_models)

    # initialize design variable array
    x0 = [copy(turbine_x_start); copy(turbine_y_start)]
    xopt = [copy(turbine_x_opt); copy(turbine_y_opt)]

    params_base = set_up_base_params(params, alpha=alpha, nrotorpoints=100)
    aep_init_base = aep_wrapper(x0, params_base)[1]
    println("aep init base: ", aep_init_base)
    aep_opt_base = aep_wrapper(xopt, params_base)[1]
    println("aep opt base: ", aep_opt_base)
    
    # get power by direction for both 1 and 100 sample points 
    turbine_power_by_direction_base_1 = turbine_power_by_direction_from_params(turbine_x_base, turbine_y_base, params)
    turbine_power_by_direction_base_100 = turbine_power_by_direction_from_params(turbine_x_base, turbine_y_base, params_base)
    turbine_power_by_direction_start_1 = turbine_power_by_direction_from_params(turbine_x_start, turbine_y_start, params)
    turbine_power_by_direction_start_100 = turbine_power_by_direction_from_params(turbine_x_start, turbine_y_start, params_base)
    turbine_power_by_direction_opt_1 = turbine_power_by_direction_from_params(turbine_x_opt, turbine_y_opt, params)
    turbine_power_by_direction_opt_100 = turbine_power_by_direction_from_params(turbine_x_opt, turbine_y_opt, params_base)

    # get wake counts
    wake_count_base = ff.wake_count_iec(turbine_x_base, turbine_y_base, wind_resource.wind_directions, rotor_diameter)
    wake_count_start = ff.wake_count_iec(turbine_x_start, turbine_y_start, wind_resource.wind_directions, rotor_diameter)
    wake_count_opt = ff.wake_count_iec(turbine_x_opt, turbine_y_opt, wind_resource.wind_directions, rotor_diameter)
    
    # save data 
    dfwcb = DataFrame(wake_count_base',:auto)
    CSV.write("turbine-wakes-ff-$case-$tuning-layout1-base.txt", dfwcb, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    dfwcs = DataFrame(wake_count_start',:auto)
    CSV.write("turbine-wakes-ff-$case-$tuning-layout$layoutid-start.txt", dfwcs, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    dfwco = DataFrame(wake_count_opt',:auto)
    CSV.write("turbine-wakes-ff-$case-$tuning-layout$layoutid-opt$n.txt", dfwco, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    dfb1 = DataFrame(turbine_power_by_direction_base_1', :auto)
    CSV.write("turbine-power-ff-1pts-$case-$tuning-layout1-base.txt", dfb1, header=string.(round.(winddirections.*180.0./pi, digits=0)))

    dfb100 = DataFrame(turbine_power_by_direction_base_100', :auto)
    CSV.write("turbine-power-ff-100pts-$case-$tuning-layout1-base.txt", dfb100, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    dfs1 = DataFrame(turbine_power_by_direction_start_1', :auto)
    CSV.write("turbine-power-ff-1pts-$case-$tuning-layout$layoutid-start.txt", dfs1, header=string.(round.(winddirections.*180.0./pi, digits=0)))

    dfs100 = DataFrame(turbine_power_by_direction_start_100', :auto)
    CSV.write("turbine-power-ff-100pts-$case-$tuning-layout$layoutid-start.txt", dfs100, header=string.(round.(winddirections.*180.0./pi, digits=0)))
    
    dfo1 = DataFrame(turbine_power_by_direction_opt_1', :auto)
    CSV.write("turbine-power-ff-1pts-$case-$tuning-layout$layoutid-opt$n.txt", dfo1, header=string.(round.(winddirections.*180.0./pi, digits=0)))

    dfo100 = DataFrame(turbine_power_by_direction_opt_100', :auto)
    CSV.write("turbine-power-ff-100pts-$case-$tuning-layout$layoutid-opt$n.txt", dfo100, header=string.(round.(winddirections.*180.0./pi, digits=0)))

end