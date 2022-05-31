using DataFrames, CSV
using PyPlot; const plt=PyPlot

function combine_datafiles(case; directory="/Users/jaredthomas/Dropbox/onedrive/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/opt-rotor-point-study/$(case)/")

    # get file names
    filenames = readdir(directory)

    # intialize dataframe
    df = DataFrame()

    # read and include all data from directory
    for file in filenames
        dfsub = DataFrame(CSV.File(directory*file, header=true))
        dfsub.layoutid = 1:length(dfsub.nrotorpoints)
        append!(df, dfsub)
    end

    # save all the results
    CSV.write(directory*"combined-results-"*case*".csv", df)

    # save results without the layouts 
    CSV.write(directory*"combined-results-no-layouts-"*case*".csv", df[!, [:id, :aepi, :aepf, :aepib, :aepfb, :aepic, :aepfc, :info, :time, :fcalls, :nrotorpoints]])

end

function load_data(case; directory="/Users/jaredthomas/Dropbox/onedrive/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/opt-rotor-point-study/$(case)/", layouts=false)
    if !layouts
        df = DataFrame(CSV.File(directory*"combined-results-no-layouts-$case.csv", header=true))
    else
        df = DataFrame(CSV.File(directory*"combined-results-$case.csv", header=true))
    end
    return df
end

function plot_data_files(case)
    #xopt,aepi,aepf,aepib,aepfb,info,out,time,fcalls,nrotorpoints
    # load data
    df = load_data(case)

    # plot average opt AEP vs bin count
    fig, ax = subplots(1)
    plt.scatter(df[!, :nrotorpoints], df[!, :aepfb])
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