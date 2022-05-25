using FLOWFarm: rotor_sample_points
using FLOWFarm; const ff=FLOWFarm 
using DelimitedFiles 
using PyPlot
using VectorizedRoutines.Matlab


# load flowfarm inputs etc 
# include("../inputfiles/model-sets/bastankhah2014-single-turb-wind-tunnel-case.jl")
# include("../inputfiles/model-sets/model_set_0_single_turbine.jl")
include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-opt.jl")

function vertical_slice(;nearwake=true)

    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_sample_points_y, rotor_sample_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(10, case="high-ti", nrotorpoints=100)
    
    wind_resource.wind_directions[1] = 3*pi/2

    # define how many points should be in the flow field
    xres = 1000
    zres = 500

    # define how far off tze ground to investigate
    maxz = 5.0*rotor_diameter[1]
    minz = 0.0
    maxx = 20.0*rotor_diameter[1]
    minx = -rotor_diameter[1]

    # set up point grid for flow field
    xrange = minx:(maxx-minx)/xres:maxx
    yrange = 0
    zrange = minz:(maxz-minz)/zres:maxz

    # run flowfarm 
    ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
    model_set, turbine_x, turbine_y, turbine_z, turbine_yaw,
    rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource)  
    
    ranges = ["res" xres zres; "max" maxx maxz; "min" minx minz]

    # save results 
    # writedlm("../../image-generation/image-data/verification/vertical-slice-ranges.txt", ranges)
    # writedlm("../../image-generation/image-data/verification/vertical-slice-interpolated.txt", ffvelocities[:, 1, :])

    # visualize 
    xg, zg = meshgrid(collect(xrange), collect(zrange))
    levels = 0:0.5:20

    plt.contourf(xg, zg, ffvelocities[:,1,:]', levels)
    plt.show()


end

function horizontal(;nearwake=true)

    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_sample_points_y, rotor_sample_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(1, case="high-ti", nrotorpoints=100)
    
    wind_resource.wind_directions[1] = 3*pi/2
    
    # define how many points should be in the flow field
    xres = 100
    yres = 100
    zres = 1
    boundary_radius = 100.0
    # define flow field domain
    maxy = boundary_radius
    miny = -boundary_radius
    maxx = boundary_radius
    minx = -boundary_radius

    # set up point grid for flow field
    xrange = minx:(maxx-minx)/xres:maxx
    yrange = miny:(maxy-miny)/yres:maxy
    zrange = hub_height[1]

    # run flowfarm 
    ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw,
        rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
        wind_resource, wind_farm_state_id=1)
    println(size(ffvelocities))
    # visualize 
    xg, yg = meshgrid(collect(xrange), collect(yrange))
    fig, ax = plt.subplots(1)
    println(size(ffvelocities))
    plt.contourf(xg, yg, ffvelocities[1,:,:])
    plt.savefig("flowfield.png")

end


 