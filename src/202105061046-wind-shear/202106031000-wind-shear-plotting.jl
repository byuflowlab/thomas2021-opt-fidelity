using FLOWFarm; const ff = FLOWFarm 
using Plots 
using DelimitedFiles

function plot_wind_shear_inflow()

    include("../inputfiles/model-sets/model_set_0_single_turbine.jl")

    # define how many points should be in the flow field
    npoints = 10

    # define how far off the ground to investigate
    maxheight = 4.0*rotor_diameter[1]./2.0

    # set up point grid for flow field
    xrange = -1*rotor_diameter[1]
    yrange = 0

    # test wind shear with uchida 2020 data for N=4 inflow
    data = readdlm("../inputfiles/results-uchida2020-shear-profiles/uchida2020-wind-shear-inflow-N4.csv", ',', skipstart=1)
    zrange = data[:,2].*rotor_diameter[1]./2.0
    inflow = data[:,1]

    hub_height = zeros(nturbines) .+ diam
    measurementheights .= diam

    shearexponent = 1.0/4.0
    wind_shear_model = ff.PowerLawWindShear(shearexponent)
    wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(100)

    ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocities = reshape(ffvelocities, (length(inflow)))

    p = scatter(inflow, zrange/(diam/2), label="Uchida 2020")
    println(size(ffvelocities))
    scatter!(p, ffvelocities./wind_speed, zrange/(diam/2), label="FLOWFarm")
    plot!(p, [1, 1], [0,4])
    plot!(p, [0.0, 1.5], [hub_height[1], hub_height[1]]./(diam/2))
    display(p)
end

function plot_wind_shear_wake_uchida()

    include("../inputfiles/model-sets/model_set_0_single_turbine.jl")

    diam = 1.0

    # set turbine design parameters
    rotor_diameter .= diam # m
    hub_height .= diam  # m
    cut_in_speed .= 3.  # m/s
    cut_out_speed .= 25.  # m/s
    rated_speed .= 11.4  # m/s
    rated_power .= 5.0E6  # W
    generator_efficiency .= 0.944

    # initialize thrust model
    ct_model = ff.ThrustModelConstantCt(0.2)
    ct_models = Vector{typeof(ct_model)}(undef, nturbines)
    for i = 1:nturbines
        ct_models[i] = ct_model
    end

    # define how many points should be in the flow field
    npoints = 10

    # define how far off the ground to investigate
    maxheight = 4.0*rotor_diameter[1]./2.0

    # set up point grid for flow field
    xrange = 30*rotor_diameter[1]
    yrange = 0

    # test wind shear with uchida 2020 data for N=4 inflow
    data = readdlm("../inputfiles/results-uchida2020-shear-profiles/uchida2020-wind-shear-wake-N4.csv", ',', skipstart=1)
    zrange = data[:,2].*rotor_diameter[1]./2.0
    inflow = data[:,1]

    hub_height .= diam
    measurementheights .= diam
    ambient_tis .= 0.00

    shearexponent = 1.0/4.0
    wind_shear_model = ff.PowerLawWindShear(shearexponent)
    wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(100)

    ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocities = reshape(ffvelocities, (length(inflow)))

    p = scatter(inflow, zrange/(diam/2), label="Uchida 2020")
    println(size(ffvelocities))
    scatter!(p, ffvelocities./wind_speed, zrange/(diam/2), label="FLOWFarm")
    plot!(p, [1, 1], [0,4])
    plot!(p, [0.0, 1.5], [hub_height[1], hub_height[1]]./(diam/2), legend=:topleft)
    display(p)
end

function plot_wind_shear_wake_bastankhah()

    include("../inputfiles/model-sets/bastankhah2014-single-turb-case-3.jl")

    # define how many points should be in the flow field
    npoints = 10

    # define how far off the ground to investigate
    maxheight = 6.0*rotor_diameter[1]./2.0

    uh = 9.0

    # set up point grid for flow field
    xrange = 7*rotor_diameter[1]
    yrange = 0

    shearexponent = 0.12539210313906432
    groundheight = 4.842460795576101
    wind_shear_model = ff.PowerLawWindShear(shearexponent, groundheight)
    wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(100)

    # test wind shear with wu 2012 data
    data = readdlm("../inputfiles/results-wu-2012-wake-profiles/wake-profile-5E-2.csv", ',', skipstart=1)
    zrange = data[:,2].*rotor_diameter[1]
    u0 = zeros(length(zrange))
    u = zeros(length(zrange))
    for i in 1:length(zrange) 
        u0[i] = ff.adjust_for_wind_shear(zrange[i], uh, 70.0, wind_shear_model.ground_height, wind_shear_model)
        u[i] = u0[i]*(1.0-data[i,1])
    end

    # test wind shear with bastankhah 2014 model data
    # data_b = readdlm("../inputfiles/results-bastankhah-2014/bastankhah-7d-5E-2.csv", ',', skipstart=1)
    data_b = readdlm("../inputfiles/results-bastankhah-2014/wake-profile-5E-2-7d-fig-7.csv", ',', skipstart=1)
    zrange_b = data_b[:,2].*rotor_diameter[1]
    u0_b = zeros(length(zrange_b))
    u_b = zeros(length(zrange_b))
    for i in 1:length(zrange_b) 
        u0_b[i] = ff.adjust_for_wind_shear(zrange_b[i], uh, 70.0, wind_shear_model.ground_height, wind_shear_model)
        u_b[i] = u0_b[i]*(1.0-data_b[i,1])
    end

    # set up test range 
    zrange_t = 0:maxheight
    u0_t = zeros(length(zrange_t))
    for i in 1:length(zrange_t) 
        u0_t[i] = ff.adjust_for_wind_shear(zrange_t[i], uh, 70.0, wind_shear_model.ground_height, wind_shear_model)
    end

    wakedeficitmodel = ff.GaussOriginal(0.04)
    model_set_bp2014 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    ffvelocitiesbp2014 = ff.calculate_flow_field(xrange, yrange, zrange_t,
        model_set_bp2014, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocitiesbp2014 = reshape(ffvelocitiesbp2014, (length(u0_t)))

    # set up wake and related models for 2016
    # wakedeficitmodel = ff.GaussYaw(0.04, 0.04, 2.32, 0.154)
    wakedeficitmodel = ff.GaussYawVariableSpread()
    model_set_2016 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    ffvelocitiesbp2016 = ff.calculate_flow_field(xrange, yrange, zrange_t,
        model_set_2016, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocitiesbp2016 = reshape(ffvelocitiesbp2016, (length(u0_t)))

    # set up wake models for Jensen Cosine
    wakedeficitmodel = ff.JensenCosine(0.075)
    
    model_set_JensenCosine = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    ffvelocitiesJensenCosine = ff.calculate_flow_field(xrange, yrange, zrange_t,
        model_set_JensenCosine, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocitiesJensenCosine = reshape(ffvelocitiesJensenCosine, (length(u0_t)))

    # set up wake models for Jensen TopHat
    wakedeficitmodel = ff.JensenTopHat(0.075)
    
    model_set_JensenTopHat = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    ffvelocitiesJensenTopHat = ff.calculate_flow_field(xrange, yrange, zrange_t,
        model_set_JensenTopHat, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocitiesJensenTopHat = reshape(ffvelocitiesJensenTopHat, (length(u0_t)))

    # set up wake models for MultiZone
    wakedeficitmodel = ff.MultiZone()
    
    model_set_JensenMultiZone = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    ffvelocitiesMultiZone = ff.calculate_flow_field(xrange, yrange, zrange_t,
        model_set_JensenMultiZone, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource)  

    ffvelocitiesMultiZone = reshape(ffvelocitiesMultiZone, (length(u0_t)))


    p = scatter((u0-u)./uh, zrange/diam, label="Wu 2012", xlabel="Normalized Wake Deficit", ylabel="z/d",legend=:topright)
    plot!(p, (u0_b-u_b)./uh, zrange_b/diam, label="Bastankhah 2014", ylim=[0, maxheight/diam])
    plot!(p, (u0_t - ffvelocitiesbp2014)./uh, zrange_t/diam, label="FLOWFarm: BP 2014")
    plot!(p, (u0_t - ffvelocitiesbp2016)./uh, zrange_t/diam, label="FLOWFarm: BP 2016")
    # plot!(p, (u0_t - ffvelocitiesJensenCosine)./uh, zrange_t/diam, label="FLOWFarm: Jensen Cosine")
    # plot!(p, (u0_t - ffvelocitiesJensenTopHat)./uh, zrange_t/diam, label="FLOWFarm: Jensen TopHat")
    # plot!(p, (u0_t - ffvelocitiesMultiZone)./uh, zrange_t/diam, label="FLOWFarm: MultiZone")
    display(p)
end