import FLOWFarm; const ff = FLOWFarm
using Statistics

function wind_farm_setup(nturbines)
    # set initial turbine x and y locations
    diam = 126.4
    data = readdlm("../inputfiles/farms/layout_38turb_round.txt",  ' ', skipstart=1)
    turbine_x = data[1:nturbines, 1].*diam
    turbine_y = data[1:nturbines, 2].*diam

    turbine_x = turbine_x .- turbine_x[1]
    turbine_y = turbine_y .- turbine_y[1]

    # set turbine base heights
    turbine_z = zeros(nturbines) .+ 0.0

    # set turbine yaw values
    turbine_yaw = zeros(nturbines)

    # set turbine design parameters
    rotor_diameter = zeros(nturbines) .+ diam # m
    hub_height = zeros(nturbines) .+ 90.0   # m
    cut_in_speed = zeros(nturbines) .+3.  # m/s
    cut_out_speed = zeros(nturbines) .+25.  # m/s
    rated_speed = zeros(nturbines) .+11.4  # m/s
    rated_power = zeros(nturbines) .+5.0E6  # W
    generator_efficiency = zeros(nturbines) .+ 1.0

    # rotor swept area sample points (normalized by rotor radius)
    nrotorpoints = 100
    # rotor_points_y, rotor_points_z = ff.rotor_sample_points(nrotorpoints, method="grid", pradius=1)
    rotor_points_y, rotor_points_z = ff.rotor_sample_points(nrotorpoints, method="sunflower", pradius=1)
    # println(size(rotor_points_y), size(rotor_points_z))
    # load tuned wind speed and ambient ti data 
    # windandtidata = readdlm("../202105181144-38-turb-tune-to-sowfa/tuned-parameters-low-ti.csv", ',', skipstart=1)
    windandtidata = readdlm("../202105181144-38-turb-tune-to-sowfa/tuned-parameters-low-ti-sowfa.csv", ',', skipstart=1)
    tunedwindspeeds = windandtidata[:, 2]

    # set flow parameters
    winddata = readdlm("../inputfiles/wind/windrose_nantucket_12dir.txt", ' ', skipstart=1)
    winddirections = winddata[:,1].*pi./180.0
    nstates = length(winddirections)
    windspeeds = zeros(nstates) .+ mean(tunedwindspeeds) #tunedwindspeeds  #zeros(nstates) .+ mean(tunedwindspeeds) #zeros(length(winddirections)) .+ 7.99 # winddata[:,2]
    windprobabilities = winddata[:,3]
    

    air_density = 1.225  # kg/m^3 (from Jen)

    tunedti = windandtidata[:,3]
    ambient_ti = mean(tunedti)
    shearexponent = 0.175 

    ambient_tis =  zeros(nstates) .+ ambient_ti #tunedti
    measurementheight = zeros(nstates) .+ 90.0 #80.0

    # load power curve
    powerdata = readdlm("../inputfiles/turbines/nrel-5mw/NREL5MWCPCT.txt", skipstart=1)
    velpoints = powerdata[:,1]
    cppoints = powerdata[:,2]

    # initialize power model
    power_model = ff.PowerModelCpPoints(velpoints, cppoints)
    power_models = Vector{typeof(power_model)}(undef, nturbines)
    for i = 1:nturbines
        power_models[i] = power_model
    end

    # load thrust curve
    ctdata = readdlm("../inputfiles/turbines/nrel-5mw/NREL5MWCPCT.txt", skipstart=1)
    velpoints = ctdata[:,1]
    ctpoints = ctdata[:,3]
    # ctpoints[ctpoints .>= 1.0] .= 0.9

    # initialize thurst model
    ct_model = ff.ThrustModelCtPoints(velpoints, ctpoints)
    # axial_induction = 1.0/3.0
    # Ct = 4.0 * axial_induction * (1.0 - axial_induction)
    # ct_model = ff.ThrustModelConstantCt(Ct)
    ct_models = Vector{typeof(ct_model)}(undef, nturbines)
    for i = 1:nturbines
        ct_models[i] = ct_model
    end

    # p = plot(velpoints, ctpoints, label="ct")
    # p = plot!(velpoints, cppoints, label="cp")
    # display(p)

    # initialize wind shear model
    wind_shear_model = ff.PowerLawWindShear(shearexponent)

    # get sorted indecies 
    sorted_turbine_index = sortperm(turbine_x)

    # initialize the wind resource definition
    wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

    # set up wake and related models
    alphastar = 2.32
    betastar = 0.154
    k1 = 0.3837
    k2 = 0.003678
    wakedeficitmodel = ff.GaussYawVariableSpread(alphastar, betastar, k1, k2, [1.0])
    # wakedeficitmodel = ff.JensenTopHat()
    wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection(alphastar, betastar, k1, k2)
    wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
    # wakecombinationmodel = ff.LinearFreestreamSuperposition()
    # wakecombinationmodel = ff.SumOfSquaresLocalVelocitySuperposition()
    localtimodel = ff.LocalTIModelMaxTI(alphastar, betastar, k1, k2)
    # localtimodel = ff.LocalTIModelNoLocalTI()

    # initialize model set
    model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    return diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set
end

diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)