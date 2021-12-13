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

    # rediscretize wind rose
    params_base.wind_resource = ff.rediscretize_windrose(params.wind_resource, 360, start=0.0, averagespeed=true)

    return params_base

end

function run_optimization(layoutid, ndirectionbins; case="high-ti", tuning="sowfa-nrel", plotresults=false, verbose=true, wec=true, nrotorpoints=1, alpha=0, savehistory=false, optimize=true, outdir="./", layoutdir="../inputfiles/farms/startinglayouts/individual/", lspacing=3.0)

    # get wind farm setup
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, _, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38, case=case, tuning=tuning, layoutid=layoutid, nrotorpoints=nrotorpoints, alpha=alpha, layoutdir=layoutdir, lspacing=lspacing, basedirs=36)

    # rediscretize wind rose
    wind_resource = ff.rediscretize_windrose(wind_resource, ndirectionbins, start=0.0, averagespeed=true)

    # get number of turbines 
    nturbines = length(turbine_x)

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
    if verbose
        println("boundary radius: $(boundary_radius)")
    end

    params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
        rotor_diameter, boundary_center, boundary_radius, obj_scale, con_scale_boundary, xyscale, hub_height, turbine_yaw, 
        ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
        wind_resource, power_models)

    # initialize design variable array
    x0 = [copy(turbine_x);copy(turbine_y)]
    
    params_base = set_up_base_params(params, alpha=alpha, nrotorpoints=100)
    aep_init_base = aep_wrapper(x0, params_base)[1]

    # print stuff if desired
    if verbose
        println("Design variables: $(size(x0)[1])")
        println("wind speed: $(mean(params.wind_resource.wind_speeds))")
        println("TI: $(mean(params.wind_resource.ambient_tis))")
        println("wind shear exp: $(params.wind_resource.wind_shear_model.shear_exponent)")
        println("AEP init (base params): $aep_init_base")
    end

    # report initial objective value
    aep_init = aep_wrapper(x0, params)[1]
    if verbose
        println("starting objective value: ", aep_init)
    end

    # add initial turbine location to plot
    if plotresults
        fig, axlayout = plt.subplots(1)
        axlayout.plot(0,0)
        for i = 1:length(turbine_x)
            axlayout.add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
        end
        # add wind farm boundary to plot
        axlayout.add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))
    
        # set up plot formatting
        axlayout.axis("square")
        axlayout.set_xlim(-boundary_radius-200,boundary_radius+200)
        axlayout.set_ylim(-boundary_radius-200,boundary_radius+200)
    end

    if !optimize 
        if plotresults
            plt.show()
        end
        return 
    end

    # set general lower and upper bounds for design variables
    lx = zeros(length(x0)) .- boundary_radius
    ux = zeros(length(x0)) .+ boundary_radius

    # set general lower and upper bounds for constraints
    ng = Int(nturbines + (nturbines)*(nturbines - 1)/2)
    lg = [-Inf*ones(Int((nturbines)*(nturbines - 1)/2)); -Inf*ones(nturbines)]
    if verbose
        println(minimum(lg), maximum(lg))
    end
    ug = [zeros(Int((nturbines)*(nturbines - 1)/2)); zeros(nturbines)]
    
    # initialize wec values and storage containers
    wec_values = [3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 1.0] 
    optaep = []
    optx = []
    opty = []
    startx = []
    starty = []
    if savehistory
        xhistory = []
    else
        xhistory = nothing
    end
    lti = nothing
    if wec
        t1 = time()
        for i = 1:length(wec_values)
            # set up models 
            if i == length(wec_values)
                model_set_local = model_set
                lti = true
                convtol = 1E-4
            else
                local_ti_model_local = ff.LocalTIModelNoLocalTI()
                model_set_local = ff.WindFarmModelSet(model_set.wake_deficit_model, model_set.wake_deflection_model, model_set.wake_combination_model, local_ti_model_local)
                model_set_local.wake_deficit_model.wec_factor[1] = wec_values[i]
                lti = false
                convtol = 1E-3
            end

            # set params 
            params = params_struct(model_set_local, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
                    rotor_diameter, boundary_center, boundary_radius, obj_scale, con_scale_boundary, xyscale, hub_height, turbine_yaw, 
                    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
                    wind_resource, power_models)

            # set SNOPT options
            snopt_opt = Dict(
                # "Derivative option" => 1,
                "Verify level" => 0,
                "Major optimality tolerance" => convtol,
                # "Major iterations limit" => 1E0,
                "Summary file" => outdir*"bins$ndirectionbins-snopt-summary-$(case)-layout-$(layoutid)-lti-$(lti)-wec-$(wec_values[i]).out",
                "Print file" => outdir*"bins$ndirectionbins-snopt-print-$(case)-layout-$(layoutid)-lti-$(lti)-wec-$(wec_values[i]).out"
            )

            # initialize solver
            solver = SNOPT(options=snopt_opt)

            # initialize SNOW options
            options = Options(;solver, derivatives=ForwardAD())

            # generate objective function 
            obj_func_wec!(g, x) = wind_farm_opt!(g, x, params, xhistory=xhistory)

            # save starting position
            push!(startx, turbine_x)
            push!(starty, turbine_y)

            # initialize starting variable array
            x0 = [copy(turbine_x);copy(turbine_y)]
            
            # optimize
            xopt, fopt, info, out = minimize(obj_func_wec!, x0, ng, lx, ux, lg, ug, options)

            # extract final turbine locations
            turbine_x = copy(xopt[1:nturbines])
            turbine_y = copy(xopt[nturbines+1:end])

            # add results to results arrays
            push!(optaep, -fopt/obj_scale)
            push!(optx, turbine_x)
            push!(opty, turbine_y)

        end
        t2 = time()
        clk = t2-t1

        # df = DataFrame(wec=wec_values, aep=optaep, xs=startx, ys=starty, x=optx, y=opty)
        # CSV.write(outdir*"optresults-wec-$case-$tuning-$layoutid.csv", df)
        # df2 = DataFrame(x=optx[end], y=opty[end])
        # CSV.write(outdir*"optresultseasy-wec-$case-$tuning-$layoutid.csv", df2)
        
    else
        xopt = nothing
        info = nothing
        out = nothing
        t1 = time()
        for i = 1:2
            if lti !== nothing
                model_set_local = model_set
                lti = true
                convtol = 1E-4
            else
                local_ti_model_local = ff.LocalTIModelNoLocalTI()
                model_set_local = ff.WindFarmModelSet(model_set.wake_deficit_model, model_set.wake_deflection_model, model_set.wake_combination_model, local_ti_model_local)
                lti = false
                convtol = 1E-3
            end

            # set params 
            params = params_struct(model_set_local, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
                    rotor_diameter, boundary_center, boundary_radius, obj_scale, con_scale_boundary, xyscale, hub_height, turbine_yaw, 
                    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
                    wind_resource, power_models)

            # set SNOPT options
            snopt_opt = Dict(
                "Derivative option" => 1,
                "Verify level" => 0,
                "Major optimality tolerance" => convtol,
                # "Major iterations limit" => 1E0,
                "Summary file" => outdir*"bins$ndirectionbins-snopt-summary-$(case)-layout-$(layoutid)-lti-$(lti).out",
                "Print file" => outdir*"bins$ndirectionbins-snopt-print-$(case)-layout-$(layoutid)-lti-$(lti).out"
            )

            # initialize solver
            solver = SNOPT(options=snopt_opt)

            # initialize SNOW options
            options = Options(;solver, derivatives=ForwardAD())

            # params.model_set.wake_deficit_model.wec_factor[1] = 1.0

            # generate wrapper function surrogate
            obj_func_no_wec!(g, x) = wind_farm_opt!(g, x, params, xhistory=xhistory)

            # optimize
            xopt, fopt, info, out = minimize(obj_func_no_wec!, x0, ng, lx, ux, lg, ug, options)
            
            # extract final turbine locations
            turbine_x = copy(xopt[1:nturbines])
            turbine_y = copy(xopt[nturbines+1:end])

            # add results to results arrays
            push!(optaep, -fopt/obj_scale)
            push!(optx, turbine_x)
            push!(opty, turbine_y)
        end
        t2 = time()
        clk = t2-t1

        # df = DataFrame(aep=optaep, xs=startx, ys=starty, x=optx, y=opty)
        # CSV.write(outdir*"optresults-$case-$tuning-$layoutid.csv", df)
        # df2 = DataFrame(x=optx[end], y=opty[end])
        # CSV.write(outdir*"optresultseasy-$case-$tuning-$layoutid.csv", df2)
    end

    # if layoutid == 1
    #     if savehistory
    #         println("calculating aep history")

    #         # initialize history plot and storage containers
    #         if plotresults
    #             fig, axhistory = plt.subplots(1)
    #         end
    #         aep_history = []
    #         xlast = []
    #         aeplast = nothing

    #         # recalculate aep for each point in the design variable history
    #         @showprogress "recalculating aep for history..." for i = 1:length(xhistory)
    #             x = xhistory[i]
    #             if xlast != x                
    #                 aeptemp = (aep_wrapper(x, params_base)/obj_scale)*1E-9
    #                 xlast = deepcopy(x)
    #                 aeplast = deepcopy(aeptemp)
    #             else
    #                 # if design variables have not changed, don't recalculate aep
    #                 aeptemp = deepcopy(aeplast)
    #             end
    #             push!(aep_history, aeptemp)
    #         end
    #         # save recalculated aep history
    #         if wec
    #             CSV.write(outdir*"aephistory-wec-$case-$tuning-$layoutid.csv", DataFrame(aep=aep_history))
    #         else
    #             CSV.write(outdir*"aephistory-$case-$tuning-$layoutid.csv", DataFrame(aep=aep_history))
    #         end
    #         # plot recalculated aep history if desired
    #         if plotresults
    #             axhistory.plot(1:length(aep_history), aep_history)
    #             axhistory.set(xlabel="Iteration", ylabel="AEP (GWh)")
    #             plt.show()
    #         end
    #     end
    # end

    # if wec && savehistory
    #     CSV.write(outdir*"loc-history-wec-$case-$tuning-$layoutid.csv", DataFrame(loc=xhistory))
    # elseif savehistory
    #     CSV.write(outdir*"loc-history-$case-$tuning-$layoutid.csv", DataFrame(loc=xhistory))
    # end

    if savehistory
        fcalls = size(xhistory)[1]
    else
        fcalls = nothing
    end
    if verbose
        println("xopt ", xopt)
    end

    aep_final = aep_wrapper(xopt, params)

    aep_final_base = aep_wrapper([copy(turbine_x);copy(turbine_y)], params_base)

    # print optimization results
    if verbose
        println("Finished in : ", clk, " (s)")
        println("info: ", info)
        # println("end objective value: ", -fopt)
        # println("major iter = ", out.major_iter)
        # println("iterations = ", out.iterations)
        # println("solve time = ", out.run_time)
        println("AEP improvement (%) = ", 100*(aep_final - aep_init)/aep_init) 
        println("Refined AEP improvement (%) = ", 100*(aep_final_base - aep_init_base)/aep_init_base)
        println("AEP init 1 pt: $aep_init")
        println("AEP final 1 pt: $aep_final")
        println("AEP init 100 pts: $aep_init_base")
        println("AEP final 100 pts: $aep_final_base")
        # println("opt locs: ", xopt)
    end

    if plotresults
        
        # add final turbine locations to plot
        for i = 1:length(turbine_x)
            axlayout.add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
        end
        
        plt.show()
    end

    return xopt, aep_init, aep_final, aep_init_base, aep_final_base, info, out, clk, fcalls
end

function run_optimization_series(nruns, case, tuning, ndirectionbins; outdir="./", layoutgen="individual", wec=true, lspacing=3.0, firstrun=1, verbose=false, plotresults=false, savehistory=true)

    layoutdir="../inputfiles/farms/startinglayouts/$layoutgen/"

    xoptdata = []
    aepinitdata = []
    aepfinaldata = []
    aepinitbasedata = []
    aepfinalbasedata = []
    infodata = []
    outdata = []
    timedata = []
    fcalldata = []

    for i = firstrun:nruns+firstrun-1
        println("running optimization $i")
        xopt, aepi, aepf, aepib, aepfb, info, out, clk, fcalls = run_optimization(i, ndirectionbins; case=case, tuning=tuning, plotresults=plotresults, verbose=verbose, wec=wec, nrotorpoints=1, alpha=0, savehistory=savehistory, optimize=true, outdir=outdir, layoutdir=layoutdir, lspacing=lspacing)
        push!(xoptdata, xopt)
        push!(aepinitdata, aepi)
        push!(aepfinaldata, aepf)
        push!(aepinitbasedata, aepib)
        push!(aepfinalbasedata, aepfb)
        push!(infodata, info)
        push!(outdata, out)
        push!(timedata, clk)
        push!(fcalldata, fcalls)
        df = DataFrame(xopt=xoptdata, aepi=aepinitdata, aepf=aepfinaldata, aepib=aepinitbasedata, aepfb=aepfinalbasedata, info=infodata, out=outdata, time=timedata, fcalls=fcalldata, ndirs=ndirectionbins)
        CSV.write(outdir*"opt-overall-results-$case-$tuning-startrun$firstrun-nruns$nruns-dirbins$ndirectionbins.csv", df)
    end

end