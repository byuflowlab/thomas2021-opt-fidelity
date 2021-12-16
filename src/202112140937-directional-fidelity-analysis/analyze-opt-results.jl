using DataFrames, CSV
using PyPlot; const plt=PyPlot

include("recalculateaep.jl")

function combine_datafiles(case; directory="/Users/jaredthomas/OneDrive - Brigham Young University/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/directional-fidelity-study/$(case)/")

    # get file names
    filenames = readdir(directory)

    # intialize dataframe
    df = DataFrame()

    # read and include all data from directory
    for file in filenames
        dfsub = DataFrame(CSV.File(directory*file, header=true))
        dfsub.layoutid = 1:length(dfsub.ndirs)
        append!(df, dfsub)
    end

    # save all the results
    CSV.write(directory*"combined-results-"*case*".csv", df)

    # save results without the layouts 
    CSV.write(directory*"combined-results-no-layouts-"*case*".csv", df[!, [:id, :aepi, :aepf, :aepib, :aepfb, :aepic, :aepfc, :info, :time, :fcalls, :ndirs]])

end

function load_data(case; directory="/Users/jaredthomas/OneDrive - Brigham Young University/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/directional-fidelity-study/$(case)/", layouts=false)
    if !layouts
        df = DataFrame(CSV.File(directory*"combined-results-no-layouts-$case.csv", header=true))
    else
        df = DataFrame(CSV.File(directory*"combined-results-$case.csv", header=true))
    end
    return df
end

function create_data_files(case)
    #xopt,aepi,aepf,aepib,aepfb,info,out,time,fcalls,ndirs
    # load data
    df = load_data(case)

    # plot average opt AEP vs bin count
    fig, ax = subplots(1)
    plt.scatter(df[!, :ndirs], df[!, :aepfb])
end

function parse_layout(dfrow)

    # get xopt original formatting
    xopt = dfrow.xopt

    # remove punctuation (except periods)
    xopt = replace.(xopt, [',','[',']']=>"")

    # split into array of strings 
    xopt = split(xopt)

    # parse
    xopt = parse.(Float64, xopt)

    # get the number of turbines
    nturbines = Int(length(xopt)/2)

    # save separately
    turbinex = xopt[1:nturbines]
    turbiney = xopt[nturbines+1:end]

    # return result
    return turbinex, turbiney
end

function recalculate_aep_nbins_100rotorsamples(;case="high-ti")

    # load existing data 
    df = load_data(case, layouts=true)

    # initialize data structure
    aepstart_nbins100rotorpoints = []
    aepopt_nbins100rotorpoints = []

    # loop through all optimizations 
    k = 0 
    @showprogress for row in eachrow(df)
        k += 1
        # parse layout x and y 
        turbine_x, turbine_y = parse_layout(row)

        # recalculate AEP 
        newaepstart, newaepend = recalculate_aep(row.layoutid, row.ndirs, turbine_x, turbine_y; case=case, tuning="sowfa-nrel", plotresults=false, verbose=true, wec=true, nrotorpoints=100, alpha=0, savehistory=false, outdir="./", layoutdir="../inputfiles/farms/startinglayouts/angle-each-circle/", lspacing=5.0)

        # save result to data structure
        push!(aepstart_nbins100rotorpoints, newaepstart)
        push!(aepopt_nbins100rotorpoints, newaepend)
        # if k==100
        #     break
        # end
    end
    plt.scatter(df.ndirs, df.aepf, label="optall-360-100")
    plt.scatter(df.ndirs, df.aepi, label="startall-360-100")
    plt.scatter(df.ndirs, df.aepfb, label="optall-ndirs-1")
    plt.scatter(df.ndirs, df.aepib, label="startall-ndirs-1")
    plt.scatter(df.ndirs[1:k], aepstart_nbins100rotorpoints, label="start ndirs 100")
    plt.scatter(df.ndirs[1:k], aepopt_nbins100rotorpoints, label="opt ndirs 100")
    plt.legend()
    # plt.show()

    # add result to dataframe 
    df.aepin100 = aepstart_nbins100rotorpoints
    df.aepfn100 = aepopt_nbins100rotorpoints

    # save dataframe to new file 
    CSV.write("combined-results-"*case*".csv", df[!, [:aepi, :aepf, :aepib, :aepfb, :info, :time, :fcalls, :ndirs, :aepin100, :aepfn100]])
    
    plt.show()
    
end