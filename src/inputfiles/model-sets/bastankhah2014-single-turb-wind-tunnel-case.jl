import FLOWFarm; const ff = FLOWFarm

turbine_x = [0.0]
turbine_y = [0.0]
turbine_z = [0.0]
turbine_yaw = [0.0]
diam = 80.0
nturbines = 1

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ diam # m
hub_height = zeros(nturbines) .+ 0.125   # m
cut_in_speed = zeros(nturbines) .+1.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+3.  # m/s
rated_power = zeros(nturbines) .+1.0E6  # W
generator_efficiency = zeros(nturbines) .+0.944

rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(100)

# initialize power model
power_model = ff.PowerModelConstantCp(0.8)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# initialize thurst model
ct_model = ff.ThrustModelConstantCt(0.42)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end

sorted_turbine_index = [1]

wind_speed = 2.2
air_density = 1.1716  # kg/m^3
# air_density = 1.2
winddirections = [270.0*pi/180.0]
windspeeds = [wind_speed]
windprobabilities = [1.0]
measurementheights = [hub_height[1]]
wtvelocities = [wind_speed]
ambient_tis = [0.07]
shearexponent = 0.08 # just a guess
wind_shear_model = ff.PowerLawWindShear(shearexponent)
wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
# wakedeficitmodel = ff.GaussOriginal(0.04)
# k1 = 
# k2 = 
# astar = 
# bstar = 

wakedeficitmodel = ff.GaussYawVariableSpread()

wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelMaxTI()

model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)