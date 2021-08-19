using ColorSchemes: length
using FLOWFarm; const ff = FLOWFarm 
using PyPlot; const plt = PyPlot
using PrettyTables
using Statistics

include("../202105111413-SOWFA-comparison/202105111413-sowfa-round-farm-comparison.jl")
include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-low-ti-alldirections.jl")

diam, turbine_x, turbine_y, turbine_z, turbine_yaw, rotor_diameter, hub_height, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, generator_efficiency, nrotorpoints, 
    rotor_points_y, rotor_points_z, winddirections, windspeeds, windprobabilities, 
    air_density, ambient_ti, shearexponent, ambient_tis, measurementheight, power_models, 
    ct_models, wind_shear_model, sorted_turbine_index, wind_resource, wakedeficitmodel, 
    wakedeflectionmodel, wakecombinationmodel, localtimodel, model_set = wind_farm_setup(38)

function load_results(case, tuning; wec=true, dir="")
    if dir==""
        dir = "opt-$case"
    end
    if wec 
        datafilename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/$dir-wec/opt-overall-results-$case-$tuning.csv"
    else
        datafilename = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/$dir-no-wec/opt-overall-results-$case-$tuning.csv"
    end
    df = DataFrame(CSV.File(datafilename, header=true))

    return df
end

function find_samples_needed(aep_samples)

    scale = 1E-0
    ns = length(aep_samples)
    aveaep = zeros(length(aep_samples))
    stdaep = zeros(length(aep_samples))
    maxaep = zeros(length(aep_samples))
    for i = 1:length(aep_samples)
        aveaep[i] = mean(aep_samples[1:i])
        stdaep[i] = std(aep_samples[1:i])
        maxaep[i] = maximum(aep_samples[1:i])
    end

    popmax = maximum(aep_samples)
    popmean = mean(aep_samples)
    popstd = std(aep_samples)

    fig, ax = subplots(1,3)

    ax[1].plot(1:ns, aveaep.*scale)
    ax[1].plot([1,ns],[0.999, 0.999].*popmean.*scale, "--k")
    ax[1].plot([1,ns],[1.001, 1.001].*popmean.*scale, "--k", label="pm 0.5%")
    ax[1].set(xlabel="Samples", ylabel="Mean Opt. AEP (GWh)")
    ax[1].legend()

    ax[2].plot(1:ns, stdaep.*scale)
    ax[2].set(xlabel="Samples", ylabel="Std. Opt. AEP (GWh)")
    ax[2].plot([1,ns],[0.95, 0.95].*popstd.*scale, "--k")
    ax[2].plot([1,ns],[1.05, 1.05].*popstd.*scale, "--k", label="pm 5%")
    ax[2].legend()

    ax[3].plot(1:ns, maxaep.*scale)
    ax[3].plot([1,ns],[0.999, 0.999].*popmax.*scale, "--k")
    ax[3].plot([1,ns],[1.001, 1.001].*popmax.*scale, "--k", label="pm 0.5%")
    ax[3].set(xlabel="Samples", ylabel="Max. Opt. AEP (GWh)")
    ax[3].legend()

    plt.legend()
    plt.tight_layout()
    plt.show()
end

function correlation_plot(df)
    dfc = df[!, isa.(eachcol(df), Vector{Float64}) .| isa.(eachcol(df), Vector{Int64})]

    nvars = length(eachcol(dfc))
    fig, ax = subplots(nvars,nvars)
    # println(size(ax))
    # println(dfc)
    # quit()
    for i = 1:2
        for j = 1:2
            println(i,j)
            if i == j 
                ax[i,j].hist(df.time) 
            else
                ax[i,j].scatter(df.fcalls, df.time)
            end
        end
    end
    
    plt.show()

end

function find_max_layout(df, case, tuning; savelayout=false)
    idx = argmax(df.aepfb)

    println("max aep from run $idx with aep of $(df.aepfb[idx])")

    if savelayout
        xopt = df.xopt[idx]
        # println(xopt)
        # remove punctuation (except periods)
        xopt = replace.(xopt, [',','[',']']=>"")
        # split into array of strings 
        xopt = split(xopt)
        # parse
        xopt = parse.(Float64, xopt)

        nturbines = Int(length(xopt)/2)

        turbinex = xopt[1:nturbines]
        turbiney = xopt[nturbines+1:end]
        
        println(typeof(turbinex))
        dfout = DataFrame(x=turbinex, y=turbiney)
        
        include("../202106021113-layouts-for-nrel/base_layouts_for_nrel.jl")

        print_optimized_layouts(df=dfout, outfile="../202106021113-layouts-for-nrel/round-38-turbines-$(case)-$(tuning)-$idx-opt.xlsx", diam=126.4)

    end
end

function calculate_results(nsamplepoints=100, case="low-ti", tuning="alldirections")

    datafile = "../../image-generation/image-data/layouts/opt/optresultsmilestone.csv"

    # extract data
    df = DataFrame(CSV.File(datafile, datarow=2, header=false))

    # name columns 
    rename!(df,:Column1 => :x,:Column2 => :y)

    turbine_powers_by_direction_ff = run_flow_farm(wind_farm_setup, x=df.x, y=df.y, use_local_ti=true, 
    nsamplepoints=nsamplepoints, alpha=0.0, verbose=false, windrose="nantucket", shearfirst=true, case=case)#ti=0.0456610699321765

    # save data 
    df = DataFrame(turbine_powers_by_direction_ff', :auto)
    CSV.write("turbine-power-ff-$(nsamplepoints)pts-$case-$tuning-opt.txt", df, header=string.(round.(winddirections.*180.0./pi, digits=0)))

end

function opt_comparison(;case="low-ti", tuning="alldirections")

    # load flowfar base data 
    basepowerfileff = "../202105111413-SOWFA-comparison/turbine-power-ff-100pts-$case-$tuning.txt"

    # load sowfa base data
    basepowerfilesowfa = "../inputfiles/results/LES/low-ti/turbine-power-low-ti.txt"

    # load flowfarm opt data 
    optpowerfileff = "turbine-power-ff-100pts-$case-$tuning-opt.txt"

    # load sowfa opt data 
    optpowerfilesowfa = "../inputfiles/results/LES/low-ti/turbine-power-low-ti-opt.txt"

    # read files to dataframes
    baseffdf = DataFrame(CSV.File(basepowerfileff, datarow=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, datarow=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, datarow=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, datarow=2, header=false))

    # compute SOWFA directional improvement
    basedirpowersowfa = sum.(eachcol(basesowfadf))
    optdirpowersowfa = sum.(eachcol(optsowfadf))
    improvementsowfa = 100.0.*(optdirpowersowfa .- basedirpowersowfa)./basedirpowersowfa

    # compute FLOWFarm directional improvement
    basedirpowerff = sum.(eachcol(baseffdf))
    optdirpowerff = sum.(eachcol(optffdf))
    improvementff = 100.0.*(optdirpowerff .- basedirpowerff)./basedirpowerff

    # compute SOFWA AEP improvement 
    aepbasesowfa = 365.0*24.0*sum(windprobabilities .* basedirpowersowfa)
    aepoptsowfa = 365.0*24.0*sum(windprobabilities .* optdirpowersowfa)

    println("SOWFA:")
    println("AEP Base: $aepbasesowfa")
    println("AEP Opt: $aepoptsowfa")
    println("AEP imp: $(100.0.*(aepoptsowfa - aepbasesowfa)/aepbasesowfa)")

    # compute FLOWFarm AEP improvement
    aepbaseff = 365.0*24.0*sum(windprobabilities .* basedirpowerff)
    aepoptff = 365.0*24.0*sum(windprobabilities .* optdirpowerff)
    
    println("FLOWFarm:")
    println("AEP Base: $aepbaseff")
    println("AEP Opt: $aepoptff")
    println("AEP imp: $(100.0.*(aepoptff - aepbaseff)/aepbaseff)")

    println("Comparison:")
    println("AEP dif base: $(100.0.*(aepbasesowfa - aepbaseff)/aepbasesowfa)")
    println("AEP dif opt: $(100.0.*(aepoptsowfa - aepoptff)/aepoptsowfa)")

    # compute differences between SOWFA and FLOWFarm
    errordf = DataFrame(BC=100.0.*(basedirpowersowfa .- basedirpowerff)./basedirpowersowfa,
                        OC=100.0.*(optdirpowersowfa .- optdirpowerff)./optdirpowersowfa)

    # compute directional difference
    dirdf = round.(DataFrame(BS=basedirpowersowfa.*1e-6, BF=basedirpowerff.*1e-6, BC=errordf.BC, 
                      OS=optdirpowersowfa.*1e-6, OF=optdirpowerff.*1e-6, OC=errordf.OC,
                      IS=improvementsowfa, IF=improvementff), digits=1)

    println(dirdf)

    table = pretty_table(dirdf; header = ["SOWFA", "BP", "Diff", "SOWFA", "BP", "Diff", "SOWFA", "BP"], backend = Val(:latex))

end
