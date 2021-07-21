using FLOWFarm; const ff=FLOWFarm
using CSV
using DataFrames
import PyPlot; const plt=PyPlot
using LsqFit

function obj_func(h, p; href=0.245, gref=0.0)
    #href=80.0, sref=8.0, hg=0.0 for sowfa
    # for Bastankhah href = 70, sref = 9.0
    # constants
    exponent = p[1]
    reference_height = href #p[3]
    reference_speed = p[2]
    ground_height = gref
    
    # get number of points for fitting
    npoints = length(h)

    # initialize shear model 
    # if p[1] < 0; p[1] = 0.0; end
    # if p[2] < 0; p[2] = 0.0; end
    shear_model = ff.PowerLawWindShear(exponent, ground_height)

    # calc v by model
    v_calc = zeros(npoints)
    for i in 1:npoints
        v_calc[i] = ff.adjust_for_wind_shear(h[i], reference_speed, reference_height, ground_height, shear_model)
    end

    return v_calc
end

function sowfa_shear(;case="high-ti")
    # set data file name for LES data
    # datafile = "../inputfiles/wind-shear-les.txt"
    
    datafile = "../202106291409-les-data-parsing/wind-shear-les-$case.txt"
    
    # load data to data frame
    df = DataFrame(CSV.File(datafile, header=0, datarow=2))
    # rename columns
    rename!(df,:Column1 => :h,:Column2 => :s)

    # optimize fit
    initial_shear = 0.1
    initial_ground = 0.0
    initial_zref = 90.0
    initial_uref = 9.0
    obj_func_local(h, p) = obj_func(h, p, href=initial_zref, gref=initial_ground)
    # fit = curve_fit(obj_func, df.h[90.0-126.4/2 .< df.h .< 90.0+126.4/2], df.s[90.0-126.4/2 .< df.h .< 90.0+126.4/2], [initial_shear, initial_ground, initial_zref, initial_uref])
    # fit = curve_fit(obj_func, df.h[5 .< df.h ], df.s[5 .< df.h ], [initial_shear, initial_ground, initial_zref, initial_uref])
    fit = curve_fit(obj_func_local, df.h[5 .< df.h ], df.s[5 .< df.h ], [initial_shear, initial_uref])
    # fit = curve_fit(obj_func, df.h, df.s, [initial_shear, initial_ground, initial_zref, initial_uref])

    # final_shear = 0.35
    final_shear = fit.param[1]
    final_ground = initial_ground
    final_zref = initial_zref
    final_uref = fit.param[2]

    # set resolution of model data
    res = 300
    # initialize output
    v_calc = zeros(res)
    # initialize heights (input)
    h_calc = collect(1:res)
    # define shear model instance
    shear_model = ff.PowerLawWindShear(final_shear, final_ground)

    # calculate speeds at each height
    for i in 1:res
        v_calc[i] = ff.adjust_for_wind_shear(h_calc[i], final_uref, final_zref, final_ground, shear_model)
    end

    # put model results in data frame
    df2 = DataFrame([h_calc, v_calc], :auto)
    rename!(df2,:x1 => :h,:x2 => :s)

    # plot LES data and model results
    fig, ax = plt.subplots()
    ax.scatter(df.s, df.h, label="LES", markersize=1)
    ax.plot([5.5,9.0],[90.0+126.4/2.0,90.0+126.4/2.0], label="Swept Area", c="Blue", linestyle="--")
    ax.plot([5.5,9.0],[90.0-126.4/2.0,90.0-126.4/2.0], c="Blue", label="", linestyle="--")
    ax.plot(df2.s, df2.h, label="Model")
    ax.set_xlabel("Wind Speed (m/s)")
    ax.set_ylabel("Height (m)")
    plt.show()

    # print optimized shear value
    println("optimized shear: ", final_shear)
    println("optimized ground: ", final_ground)
    println("optimized shear approx.: ", round(final_shear, digits=2))
    println("optimized ground approx.: ", round(final_ground, digits=2))
    println("optimized zref approx.: ", round(final_zref, digits=2))
    println("optimized uref approx.: ", round(final_uref, digits=2))

    # save model results for plotting later 
    CSV.write("wind-shear-tuned-$(case).txt", df2, header=["Height (m)", "Speed (m/s) zref=$final_zref uref=$(round(final_uref, digits=2)) shear=$(round(final_shear, digits=3)) ground=0.0"])

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

function chamorro_shear()
    # set data file name for LES data
    # datafile = "../inputfiles/wind-shear-les.txt"
    
    datafile = "../inputfiles/results-chamorro-2010/fig1windshear.csv"
    
    # load data to data frame
    df = DataFrame(CSV.File(datafile, header=0, datarow=2))
    # rename columns
    rename!(df,:Column1 => :s,:Column2 => :h)

    # optimize fit
    initial_shear = 0.9
    initial_ground = 0.0
    initial_zref = 1.0
    initial_uref = 1.0
    # fit = curve_fit(obj_func, df.h[90.0-126.4/2 .< df.h .< 90.0+126.4/2], df.s[90.0-126.4/2 .< df.h .< 90.0+126.4/2], [initial_shear, initial_ground, initial_zref, initial_uref])
    fit = curve_fit(obj_func, df.h, df.s, [initial_shear, initial_uref])
    # fit = curve_fit(obj_func, df.h, df.s, [initial_shear, initial_ground, initial_zref, initial_uref])

    # final_shear = 0.35
    final_shear = fit.param[1]
    final_ground = 0.0 #fit.param[2]
    final_zref = 0.245 #fit.param[3]
    final_uref = fit.param[2]

    # set resolution of model data
    res = 300
    # initialize output
    v_calc = zeros(res)
    # initialize heights (input)
    h_calc = collect(1/res:1/res:1)
    # define shear model instance
    shear_model = ff.PowerLawWindShear(final_shear, final_ground)

    # calculate speeds at each height
    for i in 1:res
        v_calc[i] = ff.adjust_for_wind_shear(h_calc[i], final_uref, final_zref, final_ground, shear_model)
    end

    # put model results in data frame
    df2 = DataFrame([h_calc, v_calc], :auto)
    rename!(df2,:x1 => :h,:x2 => :s)

    # plot LES data and model results
    fig, ax = plt.subplots()
    ax.scatter(df.s, df.h, label="LES")
    # ax.plot([5.5,9.0],[90.0+126.4/2.0,90.0+126.4/2.0], label="Swept Area", c="Blue", linestyle="--")
    # ax.plot([5.5,9.0],[90.0-126.4/2.0,90.0-126.4/2.0], c="Blue", label="", linestyle="--")
    ax.plot(df2.s, df2.h, label="Model")
    ax.set_xlabel("Wind Speed (m/s)")
    ax.set_ylabel("Height (m)")
    plt.show()

    # print optimized shear value
    println("optimized shear: ", final_shear)
    println("optimized ground: ", final_ground)
    println("optimized shear approx.: ", round(final_shear, digits=2))
    println("optimized ground approx.: ", round(final_ground, digits=2))
    println("optimized zref approx.: ", round(final_zref, digits=2))
    println("optimized uref approx.: ", round(final_uref, digits=2))

    # save model results for plotting later 
    # CSV.write("wind-shear-tuned-$(case).txt", df2, header=["Height (m)", "Speed (m/s) zref=$final_zref uref=$(round(final_uref, digits=2)) shear=$(round(final_shear, digits=3)) ground=0.0"])

end