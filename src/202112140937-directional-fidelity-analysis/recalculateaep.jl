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

function set_up_base_params(params; nrotorpoints=100, alpha=0, ndirectionbins=360)

    params_base = deepcopy(params)

    # set sample points 
    rotor_points_y, rotor_points_z = ff.rotor_sample_points(nrotorpoints, method="sunflower", pradius=1, alpha=alpha)

    params_base.rotor_points_y = rotor_points_y
    params_base.rotor_points_z = rotor_points_z

    # rediscretize wind rose
    params_base.wind_resource = ff.rediscretize_windrose(params.wind_resource, ndirectionbins, start=0.0, averagespeed=true)

    return params_base

end

function recalculate_aep(layoutid, ndirectionbins, turbine_xopt, turbine_yopt; case="high-ti", tuning="sowfa-nrel", plotresults=false, verbose=true, wec=true, nrotorpoints=1, alpha=0, savehistory=false, optimize=true, outdir="./", layoutdir="../inputfiles/farms/startinglayouts/angle_each_circle/", lspacing=5.0)

    # get wind farm setup
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, _, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38, case=case, tuning=tuning, layoutid=layoutid, nrotorpoints=nrotorpoints, alpha=alpha, layoutdir=layoutdir, lspacing=lspacing, basedirs=36)

    # rediscretize wind rose
    wind_resource = ff.rediscretize_windrose(wind_resource, ndirectionbins, start=0.0, averagespeed=true)

    # scale objective to be between 0 and 1
    obj_scale = 1E-9 #1E-11
    if case == "low-ti"
        obj_scale = 1E-8
    end
    con_scale_boundary = 1E-4
    xyscale = 1 #1E4

    # set wind farm boundary parameters
    boundary_center = [0.0,0.0]
    boundary_radius = 0.5 * (4000. - rotor_diameter[1])  # 1936.8 

    params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
        rotor_diameter, boundary_center, boundary_radius, obj_scale, con_scale_boundary, xyscale, hub_height, turbine_yaw, 
        ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
        wind_resource, power_models)

    params_calc = set_up_base_params(params, alpha=alpha, nrotorpoints=nrotorpoints, ndirectionbins=ndirectionbins)
    
    # initialize design variable array
    x0start = [copy(turbine_x);copy(turbine_y)]
    x0opt = [copy(turbine_xopt);copy(turbine_yopt)]
    
    aep_calc_start = aep_wrapper(x0start, params_calc)[1]
    aep_calc_opt = aep_wrapper(x0opt, params_calc)[1]

    return aep_calc_start, aep_calc_opt
end