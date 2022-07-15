using CSV: length
using SNOW
using Snopt
using DelimitedFiles 
using PyPlot
using DataFrames
using CSV
using FLOWFarm; const ff=FLOWFarm
# using Distributed
using ProgressMeter
using StatsBase
using LinearAlgebra
# import model set with wind farm and related details
include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-opt.jl")

function blockage(turbine_x, turbine_y, wind_direction, rotor_diameter, hub_height, farm_center, domain_width, domain_height; resolution=10)

    # rotate coordinates to wind direction 
    xr, yr = ff.rotate_to_wind_direction(turbine_x, turbine_y, wind_direction; center=farm_center)

    # initialize grid definition
    y = farm_center[1] - domain_width/2 : resolution : farm_center[2] + domain_width/2
    z = 0:resolution:domain_height

    # get total grid locations 
    n = length(y)*length(z)

    # initialize result grid 
    blocked_grid = zeros((length(z), length(y)))

    # loop over grid locations and sum how many are blocked 
    tiptop = hub_height + rotor_diameter/2
    for j = 1:length(z)
        if z[j] > tiptop
            continue
        end
        for i = 1:length(y)
            for k = 1:length(turbine_x)
                if norm([y[i]-yr[k], z[j]-hub_height]) <= rotor_diameter/2
                    blocked_grid[j,i] = 1
                    break
                end
            end
        end
    end

    # calculate blockage ratio 
    nblocked = sum(blocked_grid)
    blockage_ratio = nblocked/n
    println("blockage ratio: $(100*blockage_ratio)%")

    # # plot projection 
    # fig, ax = plt.subplots(1)
    # for y in yr 
    #     circle = matplotlib.patches.Circle((y, hub_height), rotor_diameter/2, color="k", fill=true)
    #     ax.add_patch(circle)
    # end

    # outflow = matplotlib.patches.Rectangle((-2500,0),5000, 1000, color="k", fill=false)
    # ax.add_patch(outflow)
    # ax.set(ylim=[0.0, 1000], xlim=[-2500, 2500], aspect="equal")
    # plt.show()

    return blockage_ratio

end

function compute_area(pos)
    x, y = (zip(pos))
    return 0.5 * abs(dot(x, roll(y, 1)) - dot(y, roll(x, 1)))
end

function blockage_all_directions(;case="high-ti", tuning="sowfa-nrel", opt=false, lspacing=3.0, alpha=0, resolution=10)
    
    # if opt
    #     layoutdir="../../image-generation/image-data/layouts/opt"
    # else
    layoutdir="../inputfiles/farms/startinglayouts/individual/"
    # end

    # get wind farm setup
    diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, _, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38, case=case, tuning=tuning, layoutid=1, nrotorpoints=nrotorpoints, alpha=alpha, layoutdir=layoutdir, lspacing=lspacing)

    # set boundary specs  
    boundary_center = [0.0, 0.0]
    domain_width = 5000
    domain_height = 1000
    les_radius = 2500.0
    boundary_radius = les_radius - 500 - rotor_diameter[1]/2

    println("Base case")
    # initialize data containers
    blockage_ratio_base = zeros(length(wind_resource.wind_directions))

    # get blockage values for all direction 
    for i = 1:length(wind_resource.wind_directions)
        blockage_ratio_base[i] = blockage(turbine_x, turbine_y, wind_resource.wind_directions[i], rotor_diameter[1], hub_height[1], boundary_center, domain_width, domain_height, resolution=resolution)
    end

    if case == "high-ti"
        n = 4
        layoutid = 83
    elseif case == "low-ti"
        n = 3
        layoutid = 385
    end
    optdatafile = "../../image-generation/image-data/layouts/opt/$case-$tuning-opt$n-layout$layoutid-aec-wec.csv"
    data = DataFrame(CSV.File(optdatafile))#.*rotor_diameter[1] .+ (2500 - boundary_radius)
    turbine_x = data[:,1]
    turbine_y = data[:,2]
    
    println("Optimized case")

    # initialize data containers
    blockage_ratio_opt = zeros(length(wind_resource.wind_directions))

    # get blockage values for all direction 
    for i = 1:length(wind_resource.wind_directions)
        blockage_ratio_opt[i] = blockage(turbine_x, turbine_y, wind_resource.wind_directions[i], rotor_diameter[1], hub_height[1], boundary_center, domain_width, domain_height, resolution=resolution)
    end

    # save results 
    # save results 
    df = DataFrame(direction=wind_resource.wind_directions, blockage=blockage_ratio_base)
    CSV.write("base_case_blockage_$(case)_base_res$resolution.csv", df)
    df = DataFrame(direction=wind_resource.wind_directions, blockage=blockage_ratio_opt)
    CSV.write("base_case_blockage_$(case)_opt_res$resolution.csv", df)

    # test
    test_ratio = (pi*(126.4/2)^2)/(5000*1000)
    turbine_x = [0.0]
    turbine_y = [0.0]

    # get blockage values for test
    blockage_ratio_test = blockage(turbine_x, turbine_y, wind_resource.wind_directions[1], rotor_diameter[1], hub_height[1], boundary_center, domain_width, domain_height, resolution=resolution)
   

    println("Test Ratio Analytic: $test_ratio")
    println("Test Ratio Numeric: $blockage_ratio_test)")
    println("Error: $(test_ratio .- blockage_ratio_test)")

end