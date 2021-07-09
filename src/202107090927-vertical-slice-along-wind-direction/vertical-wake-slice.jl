using FLOWFarm: rotor_sample_points
using FLOWFarm; const ff=FLOWFarm 
using DelimitedFiles 
using PyPlot
using VectorizedRoutines.Matlab


# load flowfarm inputs etc 
include("../inputfiles/model-sets/bastankhah2014-single-turb-wind-tunnel-case.jl")

function vertical_slice(;nearwake=true)
    
    # define how many points should be in the flow field
    xres = 1000
    zres = 500

    # define how far off tze ground to investigate
    maxz = 2.0*rotor_diameter[1]
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
    rotor_diameter, hub_height, sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource, shearfirst=false)   
    
    ranges = ["res" xres zres; "max" maxx maxz; "min" minx minz]

    # save results 
    writedlm("../../image-generation/image-data/verification/vertical-slice-ranges.txt", ranges)
    writedlm("../../image-generation/image-data/verification/vertical-slice-interpolated.txt", ffvelocities[:, 1, :])

    # visualize 
    # xg, zg = meshgrid(collect(xrange), collect(zrange))
    # levels = 1.5:0.1:2.4
    # plt.contourf(xg, zg, ffvelocities[:,1,:], levels)
    # plt.show()


end


 