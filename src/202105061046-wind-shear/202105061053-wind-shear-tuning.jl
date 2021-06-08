using FLOWFarm; const ff=FLOWFarm
using CSV
using DataFrames
using Plots
using StatsPlots
using LsqFit

function obj_func(h, p; href=70.0, sref=9.0)
    #href=80.0, sref=8.0, hg=0.0 for sowfa
    # constants
    reference_height = href
    reference_speed = sref
    ground_height = p[2]

    # get number of points for fitting
    npoints = length(h)

    # initialize shear model 
    if p[1] < 0; p = 0.0; end
    shear_model = ff.PowerLawWindShear(p[1])

    # calc v by model
    v_calc = zeros(npoints)
    for i in 1:npoints
        v_calc[i] = ff.adjust_for_wind_shear(h[i], reference_speed, reference_height, ground_height, shear_model)
    end
    return v_calc
end

function sowfa_shear()
    # set data file name for LES data
    datafile = "../inputfiles/wind-shear-les.txt"
    # load data to data frame
    df = DataFrame(CSV.File(datafile, header=0, datarow=2))
    # rename columns
    rename!(df,:Column1 => :h,:Column2 => :s)

    # optimize fit
    initial_shear = 0.1
    fit = curve_fit(obj_func, df.h, df.s, [initial_shear])

    # final_shear = 0.35
    final_shear = fit.param[1]

    # set resolution of model data
    res = 250
    # initialize output
    v_calc = zeros(res)
    # initialize heights (input)
    h_calc = collect(1:res)
    # define shear model instance
    shear_model = ff.PowerLawWindShear(final_shear)

    # calculate speeds at each height
    for i in 1:res
        v_calc[i] = ff.adjust_for_wind_shear(h_calc[i], df.s[4], df.h[4], 0.0, shear_model)
    end

    # put model results in data frame
    df2 = DataFrame([h_calc, v_calc], :auto)
    rename!(df2,:x1 => :h,:x2 => :s)

    # plot LES data and model results
    p = scatter(df.s, df.h, label="LES")
    plot!(df2.s, df2.h, label="Model")

    # print optimized shear value
    println("optimized shear: ", final_shear)
    println("optimized shear approx.: ", round(final_shear, digits=2))
end

function bastankhah_shear()

    # optimized shear: 0.12539210313906432
    # optimized shear approx.: 0.13
    # optimized ground height: 4.842460795576101
    # set data file name for LES data
    datafile = "../inputfiles/results-wu-2012-wake-profiles/shear-profile-5E-2.csv"
    # datafile = "../inputfiles/porteagel2011/porteagel2011-shear-fig-12a.csv"
    # load data to data frame
    df = DataFrame(CSV.File(datafile, header=0, datarow=2))
    # rename columns
    rename!(df,:Column1 => :s,:Column2 => :h)
    df.h .*= 80.0
    # optimize fit
    initial_shear = 0.0

    fit = curve_fit(obj_func, df.h[df.h .> 10], df.s[df.h .> 10], [initial_shear, 0.0])

    # final_shear = 0.35
    final_shear = fit.param[1]
    ground_height = fit.param[2]

    # set resolution of model data
    res = 250
    # initialize output
    v_calc = zeros(res)
    # initialize heights (input)
    h_calc = collect(175/res:175/res:175)
    # define shear model instance
    shear_model = ff.PowerLawWindShear(final_shear, ground_height)

    # calculate speeds at each height
    for i in 1:res
        v_calc[i] = ff.adjust_for_wind_shear(h_calc[i], 9.0, 70.0, ground_height, shear_model)
    end

    # put model results in data frame
    df2 = DataFrame([h_calc, v_calc], :auto)
    rename!(df2,:x1 => :h,:x2 => :s)

    # plot LES data and model results
    p = scatter(df.s, df.h, label="LES")
    plot!(df2.s, df2.h, label="Model", legend=:topleft)
    display(p)

    # print optimized shear value
    println("optimized shear: ", final_shear)
    println("optimized shear approx.: ", round(final_shear, digits=2))
    println("optimized ground height: $ground_height")
end