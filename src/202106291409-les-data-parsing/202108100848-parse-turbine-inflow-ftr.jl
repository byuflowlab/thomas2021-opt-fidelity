using Arrow
using DelimitedFiles 
using DataFrames
using CSV
using Statistics
using Interpolations
using FLOWFarm; const ff = FLOWFarm

function read_les_data(case; dframe=true, save=false)
    if case == "low-ti"
        # filename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/low_TI/lowTI_base_turbAvgVels.ftr"
        filename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/low_TI/lowTI_base_turbAvgVels_Wind1Velx.ftr"
    else
        # filename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/high_TI/highTI_base_turbAvgVels.ftr"
        filename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/high_TI/highTI_base_turbAvgVels_Wind1Velx.ftr"
    end
    
    table = Arrow.Table(filename)

    dfraw = DataFrame(table)

    re = r"(?<=wd)\d+"

    # initialize container array for data 
    data = zeros((38,12))
    winddirections = zeros(12)
    # parse data into array
    for i = 1:length(dfraw.case)
        st = dfraw.case[i]
        winddirections[i] = parse(Int64, match(re,st).match)
        data[:, i] .= dfraw.turb_u_avg[i]
    end

    println(winddirections)
    df = DataFrame(data, :auto)
    rename!(df, string.(winddirections))

    if save
        CSV.write("wind-shear-les-$case.txt", df)
    end

    if dframe
        return df
    else
        return data, winddirections
    end
end

function get_inflow_u_ave(df::DataFrame, case="high-ti")
    df_plane = df[(df[:, :x].==5.0),:]
    h = unique(df_plane[(df_plane[:, :z].<300.0), :z])
    u = zeros(0)
    for val in h
        tmp_avg = mean(df_plane[(df_plane[:, :z].==val),:u])
        println(val, tmp_avg)
        push!(u, tmp_avg)
    end
    
    df_out = DataFrame(h=h, u=u)

    CSV.write("wind-shear-les-$case.txt", df_out)
end

function get_turbine_inflow_u_ave(case="high-ti")

    # set turbine specs 
    rotordiameter = 126.4 
    hubheight = 90.0

    # load wind rose information
    wind_rose = readdlm("../inputfiles/wind/windrose_nantucket_12dir.txt", skipstart=1)

    # set wind directions
    winddirections = wind_rose[:,1]

    # load base wind turbine locations
    layout = readdlm("../inputfiles/farms/layout_38turb_round.txt", skipstart=1)

    # adjust locations to appropriate wind turbine diameter 
    layout .*= rotordiameter

    # extract layout x and y
    turbinex = layout[:, 1]
    turbiney = layout[:, 2]

    # shift turbine locations to center middle turbine at (2500, 2500) as in LES
    turbinex .= .- turbinex[1] .+ 2500.0
    turbiney .= .- turbiney[1] .+ 2500.0

    # set up storage for results 
    uave = zeros((12,38))

    # loop over all wind directions
    for i = 1:12

        wd = winddirections[i]

        # rotate turbines to current wind direction
        xtemp, ytemp = ff.rotate_to_wind_direction(turbinex, turbiney, wd)

        # load flow field for current wind direction
        flowfield = read_les_data(case, wd=Int(wd), dframe=true)

        # set up interpolation object on u speed
        itp = interpolate(df[:, 1:3], df[:, 4])

        # set up sample locations (sunflower) for a single turbine

        # run interpolation for u-ave at each turbine in each direction (loop)
    
    end
    # CSV.write("turbine-inflow-u-ave-by-direction-$case.txt", df_out)
end

function parse_ftr(case)
    df = read_les_data(case)
    println(df)
    println(maximum(maximum(eachcol(df))))
    # get_inflow_u_ave(df, case)
end

