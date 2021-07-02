using SNOW
using Snopt
using DelimitedFiles 
using PyPlot
# using Distrib@uted

# import model set with wind farm and related details
include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-low-ti.jl")

# set globals struct for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_center
    boundary_radius
    obj_scale
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

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    boundary_center = params.boundary_center
    boundary_radius = params.boundary_radius

    # find the number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundary_center, boundary_radius, turbine_x, turbine_y)
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
function wind_farm_opt!(g, x, params)

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x, params)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x, params)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[1:(end-nturbines)] = spacing_con[:]
    g[end-nturbines+1:end] = boundary_con[:]
    
    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x, params)[1]
    
    return AEP #, dAEP_dx, dcdx, fail
end

function run_single_optimization(layout; plotresults=true, verbose=true)

    # get wind farm setup
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)

    # get number of turbines 
    nturbines = length(turbine_x)
    # scale objective to be between 0 and 1
    obj_scale = 1E-7#1E-11
    xyscale = 1 #1E4

    # set wind farm boundary parameters
    boundary_center = [0.0,0.0]
    boundary_radius = 0.5 * (4000. - rotor_diameter[1])  # 1936.8 
    if verbose
        println("boundary radius: $(boundary_radius)")
    end

    params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
        rotor_diameter, boundary_center, boundary_radius, obj_scale, xyscale, hub_height, turbine_yaw, 
        ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
        wind_resource, power_models)

    # initialize design variable array
    x0 = [copy(turbine_x);copy(turbine_y)]
    xinit = deepcopy(x0)
    if verbose
        println(size(x0))
    end
    # report initial objective value
    aep_init = aep_wrapper(x0, params)[1]
    if verbose
        println("starting objective value: ", aep_init)
    end

    if plotresults
        plot(0,0)
        # add initial turbine location to plot
        for i = 1:length(turbine_x)
            plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
        end
    end

    # generate wrapper function surrogate
    obj_func!(g, x) = wind_farm_opt!(g, x, params)

    # set SNOPT options
    snopt_opt = Dict(
        "Derivative option" => 1,
        "Verify level" => 1,
        "Major optimality tolerance" => 1e-6,
        # "Major iterations limit" => 1E0,
        "Summary file" => "snopt_summary.out",
        "Print file" => "snopt_print.out"
    )

    # initialize solver
    solver = SNOPT(options=snopt_opt)

    # initialize SNOW options
    options = Options(;solver, derivatives=ForwardAD())

    # set general lower and upper bounds for design variables
    lx = zeros(length(x0)) .- boundary_radius
    ux = zeros(length(x0)) .+ boundary_radius

    # set general lower and upper bounds for constraints
    ng = Int(nturbines + (nturbines)*(nturbines - 1)/2)
    # println("ng: ", ng)
    lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
    if verbose
        println(minimum(lg), maximum(lg))
    end
    # quit()
    ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]
    # println("ug: ", minimum(lg), maximum(ug))
    # run and time optimization
    # quit()
    t1 = time()

    # println(length(x0))
    # println(length(lg))
    xopt, fopt, info, out = minimize(obj_func!, x0, ng, lx, ux, lg, ug, options)
    t2 = time()
    clk = t2-t1

    if verbose
        println("xopt ", xopt)
    end
    aep_final = aep_wrapper(xopt, params)

    # print optimization results
    if verbose
        println("Finished in : ", clk, " (s)")
        println("info: ", info)
        println("end objective value: ", -fopt)
        # println("major iter = ", out.major_iter)
        # println("iterations = ", out.iterations)
        # println("solve time = ", out.run_time)
        println("AEP improvement (%) = ", 100*(aep_final - aep_init)/aep_init) 
        println("opt locs: ", xopt)
    end
    # extract final turbine locations
    turbine_x = copy(xopt[1:nturbines])
    turbine_y = copy(xopt[nturbines+1:end])

    if plotresults
        # add final turbine locations to plot
        for i = 1:length(turbine_x)
            plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
        end
        
        # add wind farm boundary to plot
        plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))
    
        # set up and show plot
        axis("square")
        xlim(-boundary_radius-200,boundary_radius+200)
        ylim(-boundary_radius-200,boundary_radius+200)
        plt.show()
    end

    return xopt, fopt, info, out
end