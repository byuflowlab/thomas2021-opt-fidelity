using FLOWFarm: wake_combination_model
using Base: NonReshapedReinterpretArray, front
import PyPlot; const plt = PyPlot
using DataFrames 
using CSV
using Colors, ColorSchemes
using FLOWFarm; const ff = FLOWFarm
using PrettyTables
using DelimitedFiles
using LaTeXStrings
import VectorizedRoutines; const ml = VectorizedRoutines.Matlab

function custum_color_map(;idx=[3, 1, 4])
    colors = [colorant"#BDB8AD", colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    # @pyimport matplotlib.colors as matcolors
    # cmap = matcolors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],"A")

    return plt.ColorMap("BlueGrayOrange", colors[idx])
end

function heatmap(data, row_labels, col_labels; ax=nothing, cbar_kw=Dict(), cbarlabel="", use_cbar=true, labelpixels=true, vcolor="w", edgecolor="w", boxmaxmin=true, edgecolors=nothing)
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the pcolormesh call.
    """

    if ax === nothing
        ax = plt.gca()
    end

    # find max and min cells 
    maxloc = findmax(data)[2]
    minloc = findmin(data)[2]

    # label the pixels
    if edgecolors === nothing
        edgecolors = []
        if boxmaxmin
            for i = 1:length(row_labels)
                for j = 1:length(col_labels)
                    if (i == maxloc[1] && j == maxloc[2])
                        push!(edgecolors, "k")
                    elseif (i == minloc[1] && j == minloc[2])
                        push!(edgecolors, "k")
                    else
                        push!(edgecolors, edgecolor)
                    end
                end
            end
        end
    end

    # Plot the heatmap
    vmax=round(maximum(abs.(cbar_kw[:ticks])), digits=0)
    vmin = -vmax
    
    im = ax.pcolormesh(data, edgecolor=edgecolors, cmap=cbar_kw[:cmap], vmin=vmin, vmax=vmax)

    # Create colorbar
    if use_cbar
        cbar = ax.figure.colorbar(im; ax=ax, cbar_kw...)
        cbar.ax.set_xlabel(cbarlabel, rotation=0, va="bottom", labelpad=16)
    else
        cbar = nothing
    end

    # make pixels square
    ax.set(aspect="equal")

    # set tick labels and locations for x axis
    ax.set(xticklabels=col_labels, xticks=(1:length(col_labels)).-0.5)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=true, bottom=false, labeltop=true, labelbottom=false)

    # set tick labels and locations for y axis
    ax.set(yticklabels=row_labels, yticks=(1:length(row_labels)).-0.5)

    # reverse y axis to make the result more table-like
    ax.invert_yaxis()

    # remove ticks
    ax.tick_params(which="minor", top=false, right=false, bottom=false, left=false)
    ax.tick_params(which="major", top=false, right=false, bottom=false, left=false)

    # removes spines
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.spines["left"].set_visible(false)

    # label the pixels
    if labelpixels
        for i = 1:length(row_labels)
            for j = 1:length(col_labels)
                ax.text(j-0.5,i-0.5,data[i,j],
                        ha="center",va="center",color=vcolor)
            end
        end
    end

    return im, cbar
end

function wind_shear_tuning(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="windshear-", case="high-ti")

    # load data 
    df_model = DataFrame(CSV.File("./image-data/wind-shear/$(case)/wind-shear-tuned-$(case).txt", datarow=2, header=false))
    df_les = DataFrame(CSV.File("./image-data/wind-shear/$(case)/wind-shear-les-$(case).txt", datarow=2, header=false))

    # name columns 
    rename!(df_model,:Column1 => :h,:Column2 => :s)
    rename!(df_les,:Column1 => :h,:Column2 => :s)

    
    # set variables 
    # zref=90.0 
    # uref=7.87
    # shear=0.091
    # ground=0.0

    maxv = 10#round(maximum(df_model.s), digits=0)
    minv = 4#round(minimum(df_model.s), digits=0)
    rotor_diameter = 126.4 # m 
    hub_height = 90.0 # m 
    swept_top = hub_height + rotor_diameter/2.0
    swept_bottom = hub_height - rotor_diameter/2.0

    # create plot 
    fig, ax = plt.subplots(figsize=(5.5,2.5))

    # plot les data 
    ax.scatter(df_les.s, df_les.h, color=colors[2], s=10)

    if case=="high-ti"
        ax.annotate("LES", (8.9, 130), color=colors[2], alpha=1.0, size=fontsize)
    elseif case=="low-ti"
        ax.annotate("LES", (8.4, 130), color=colors[2], alpha=1.0, size=fontsize)
        # ax.annotate("LES", (7.9, 125), color=colors[2], alpha=1.0, size=fontsize)
    end

    # plot model 
    ax.plot(df_model.s, df_model.h, color=colors[3])

    if case=="high-ti"
        ax.annotate("Curve Fit", (7.5, 130), color=colors[3], alpha=1.0, size=fontsize)
    elseif case=="low-ti"
        ax.annotate("Curve Fit", (7.5, 130), color=colors[3], alpha=1.0, size=fontsize)
        # ax.annotate("Curve Fit", (8.25, 95), color=colors[3], alpha=1.0, size=fontsize) 
    end
    # add rotor swept area 
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_top, swept_top], color=colors[1], linestyle="--")
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_bottom, swept_bottom], color=colors[1], linestyle="--")
    ax.fill_between(minv:maxv, swept_bottom,swept_top, alpha=0.15, color=colors[1])
    ax.annotate("Rotor-Swept Region", (5.5, (hub_height+swept_bottom)/2.0), color=colors[1], alpha=1.0, size=fontsize)

    # add hub height
    plt.plot([minv,maxv], [hub_height, hub_height], color=colors[1], linestyle="-.")

    # add axis labels
    plt.xlabel("Wind Speed (m/s)", fontsize=fontsize)
    plt.ylabel("Height (m)", fontsize=fontsize)

    # remove top and right spines 
    # Hide the right and top spines
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    # set limits 
    ax.set_xlim([minv, maxv])
    ax.set_ylim([0.0, 310.0])

    # set ticks 
    plt.xticks(minv:maxv, size=fontsize)
    plt.yticks(0:100:300, size=fontsize)

    # make sure everything fits 
    plt.tight_layout()

    # save figure
    if showfigs
        show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*case*".pdf", transparent=true)
    end
end

function layout(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="layout-", case="low-ti-opt")
    
    # load data
    if case == "low-ti-opt"
        datafile = "image-data/layouts/opt/optresultsmilestone.csv"
    else
        ErrorException("Case not available")
    end

    # extract data
    df = DataFrame(CSV.File(datafile, datarow=2, header=false))

    # name columns 
    rename!(df,:Column1 => :x,:Column2 => :y)

    # set up plot 
    fig, ax = plt.subplots() 

    diam = 126.4
    les_radius = 2500.0
    farm_radius = les_radius - 500 - diam/2

    # plot using FLOWFarm
    ff.plotlayout!(ax, df.x .+ les_radius, df.y .+ les_radius, zeros(length(df.x)).+diam)

    # add boundary 
    circle = matplotlib.patches.Circle((les_radius, les_radius), farm_radius, fill=nothing, color=colors[1])
    ax.add_patch(circle)

    # remove top and right spines 
    # Hide the right and top spines
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["left"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.axis("off")

    # set limits 
    ax.set_xlim([0.0, les_radius*2.0])
    ax.set_ylim([0.0, les_radius*2.0])

    # adjust aspeck to equal 
    ax.set(aspect="equal")

    # make sure everything fits 
    plt.tight_layout()

    
    # save figure
    if showfigs
        show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*case*".pdf", transparent=true)
    end
end

function opt_comparison_table()

    # flowfarm base data 
    basepowerfileff = "image-data/power/turbine-power-low-ti-ff-100pts.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/turbine-power-low-ti-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/turbine-power-low-ti-ff-100pts-opt.txt"

    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/turbine-power-low-ti-les-opt.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    baseffdf = DataFrame(CSV.File(basepowerfileff, datarow=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, datarow=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, datarow=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, datarow=2, header=false))
    winddf = DataFrame(CSV.File(winddatafile, datarow=2, header=false))

    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    # compute SOWFA directional improvement
    basedirpowersowfa = sum.(eachcol(basesowfadf))
    optdirpowersowfa = sum.(eachcol(optsowfadf))
    improvementsowfa = 100.0.*(optdirpowersowfa .- basedirpowersowfa)./basedirpowersowfa

    # compute FLOWFarm directional improvement
    basedirpowerff = sum.(eachcol(baseffdf))
    optdirpowerff = sum.(eachcol(optffdf))
    improvementff = 100.0.*(optdirpowerff .- basedirpowerff)./basedirpowerff

    # compute SOFWA AEP improvement 
    aepbasesowfa = 365.0*24.0*sum(winddf.p .* basedirpowersowfa)
    aepoptsowfa = 365.0*24.0*sum(winddf.p .* optdirpowersowfa)

    println("SOWFA:")
    println("AEP Base: $aepbasesowfa")
    println("AEP Opt: $aepoptsowfa")
    aepimpsowfa = 100.0.*(aepoptsowfa - aepbasesowfa)/aepbasesowfa
    println("AEP imp: $(aepimpsowfa)")

    # compute FLOWFarm AEP improvement
    aepbaseff = 365.0*24.0*sum(winddf.p .* basedirpowerff)
    aepoptff = 365.0*24.0*sum(winddf.p .* optdirpowerff)
    
    println("FLOWFarm:")
    println("AEP Base: $aepbaseff")
    println("AEP Opt: $aepoptff")
    aepimpff = 100.0.*(aepoptff - aepbaseff)/aepbaseff
    println("AEP imp: $(aepimpff)")

    println("Comparison:")
    aepdiffbase = 100.0.*(aepbasesowfa - aepbaseff)/aepbasesowfa
    aepdiffopt = 100.0.*(aepoptsowfa - aepoptff)/aepoptsowfa
    println("AEP dif base: $(aepdiffbase)")
    println("AEP dif opt: $(aepdiffopt)")

    # compute differences between SOWFA and FLOWFarm
    errordf = DataFrame(BC=100.0.*(basedirpowersowfa .- basedirpowerff)./basedirpowersowfa,
                        OC=100.0.*(optdirpowersowfa .- optdirpowerff)./optdirpowersowfa)

    # scale and round data 
    basedirpowersowfa = round.(basedirpowersowfa.*1e-6, digits=1)
    basedirpowerff = round.(basedirpowerff.*1e-6, digits=1)
    optdirpowersowfa = round.(optdirpowersowfa.*1e-6, digits=1)
    optdirpowerff = round.(optdirpowerff.*1e-6, digits=1)

    # add units
    dirdata = []
    for i = 1:length(winddf.d)
        bs = basedirpowersowfa[i]
        bf = basedirpowerff[i]
        bc = round(errordf.BC[i],digits=1)
        os = optdirpowersowfa[i]
        of = optdirpowerff[i]
        oc = round(errordf.OC[i], digits=1)
        is = round(improvementsowfa[i], digits=1)
        iff = round(improvementff[i], digits=1)
        
        println("\\multicolumn{1}{l}{\$$(Int(round(winddf.d[i], digits=0)))^{\\circ}\$} & &")
        println("\\SI[per-mode=symbol]{$bs}{\\mega\\watt} &")
        println("\\SI[per-mode=symbol]{$bf}{\\mega\\watt} &")
        println("\\SI[per-mode=symbol]{$bc}{\\percent} & &")
        println("\\SI[per-mode=symbol]{$os}{\\mega\\watt} &")
        println("\\SI[per-mode=symbol]{$of}{\\mega\\watt} &")
        println("\\SI[per-mode=symbol]{$oc}{\\percent} & &")
        println("\\SI[per-mode=symbol]{$is}{\\percent} &")
        println("\\SI[per-mode=symbol]{$iff}{\\percent} \\\\")
        
        # println(dirdata[i,:])
    end

    println("\\multicolumn{1}{l}{AEP} & &")
    println("\\SI[per-mode=symbol]{$(round(aepbasesowfa*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    println("\\SI[per-mode=symbol]{$(round(aepbaseff*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    println("\\SI[per-mode=symbol]{$(round(aepdiffbase, digits=1))}{\\percent} & &")
    println("\\SI[per-mode=symbol]{$(round(aepoptsowfa*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    println("\\SI[per-mode=symbol]{$(round(aepoptff*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    println("\\SI[per-mode=symbol]{$(round(aepdiffopt, digits=1))}{\\percent} & &")
    println("\\SI[per-mode=symbol]{$(round(aepimpsowfa, digits=1))}{\\percent} &")
    println("\\SI[per-mode=symbol]{$(round(aepimpff, digits=1))}{\\percent} \\\\")

end

function directional_comparison_figure(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="directional-comparison")
    
    # flowfarm base data 
    basepowerfileff = "image-data/power/turbine-power-ff-100pts-low-ti-alldirections.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/turbine-power-low-ti-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/turbine-power-ff-100pts-low-ti-alldirections-opt.txt"

    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/turbine-power-low-ti-les-opt.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    baseffdf = DataFrame(CSV.File(basepowerfileff, datarow=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, datarow=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, datarow=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, datarow=2, header=false))
    winddf = DataFrame(CSV.File(winddatafile, datarow=2, header=false))

    println("compare")
    println(sum.(eachcol(baseffdf .- basesowfadf)))
    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    # compute directional data 
    basedirpowerff = sum.(eachcol(baseffdf))
    optdirpowerff = sum.(eachcol(optffdf))
    basedirpowersowfa = sum.(eachcol(basesowfadf))
    optdirpowersowfa = sum.(eachcol(optsowfadf))

    # create directional power bar charts
    fig, ax = plt.subplots(1, figsize=[6,4])

    # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
    ax.plot(winddf.d , optdirpowersowfa.*1e-6, color=colors[3], marker="o", linestyle="--")
    ax.annotate("SOWFA-Optimized", (160, 65), color=colors[3])
    # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
    ax.plot(winddf.d , basedirpowersowfa.*1e-6, color=colors[2], marker="o", linestyle="--")
    ax.annotate("SOWFA-Base", (160, 55), color=colors[2])
    # ax.bar(winddf.d .+ 5, optdirpowerff.*1E-6, label="BP-Opt", width=10, color=colors[4])
    ax.plot(winddf.d , optdirpowerff.*1e-6, color=colors[4], marker="o", linestyle="--")
    ax.annotate("BP-Optimized", (160, 59), color=colors[4])
    # ax.bar(winddf.d .+ 5, basedirpowerff.*1E-6, label="BP-Base", width=10, color=colors[1])
    ax.plot(winddf.d , basedirpowerff.*1e-6, color=colors[1], marker="o", linestyle="--")
    ax.annotate("BP-Base", (160, 50), color=colors[1])

    # format the figure
    ax.set(xticks=winddf.d, ylim=[40, 70], xlabel="Direction (deg.)", ylabel="Directional Power (MW)")
    ax.legend(frameon=false,ncol=2)

    # remove upper and right bounding box
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    # make everything fit
    plt.tight_layout()

    # save figure
    if showfigs
        show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*".pdf", transparent=true)
    end
end

function turbine_comparison_figures(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="turbine-comparison", case="low-ti", tuning="alldirections")
    
    function plot_turbine_heatmap(data, winddirections, vmin, vmax; wake_count=nothing)

        # number of turbines in the farm
        nturbines = 38

        # number of wind directions 
        ndirections = length(winddirections)

        # intialize figure
        fig, ax = plt.subplots(figsize=(8,4))

        # set tick locations and labels
        ticks = vmin:4:vmax

        # determine edge colors based on wake count 
        if wake_count === nothing
            edgecolors = nothing
        else
            # set colors
            front_color = "w"
            back_color = colors[3]
            mid_color = colors[1]
            # initialize color array 
            edgecolors = []
            for i = 1:ndirections 
                for j = 1:nturbines 
                    if wake_count[i,j] == 0 
                        push!(edgecolors, front_color)
                    elseif wake_count[i,j] > 2
                        push!(edgecolors, back_color)
                    else
                        push!(edgecolors, mid_color)
                    end
                end
            end
        end

        # generate a custom color map 
        cmap = custum_color_map()

        # generate dict of colormap options
        d = Dict(:shrink => 0.47, :ticks=>ticks, :aspect=>20, :orientation=>"horizontal", :cmap=>cmap, :pad=>0.05)
        
        # set row labels
        rowlabels = convert.(Int64, round.((winddirections), digits=0))

        # create heatmap
        im, cbar = heatmap(data, rowlabels, 1:nturbines, ax=ax, edgecolors = edgecolors,
                cbarlabel="Turbine Power Error as Percent of Max SOWFA Turbine Power", cbar_kw=d)

        # remove upper and right bounding box
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        # make everything fit
        plt.tight_layout()

        return ax
    end

    function errors(first, second; method="absolute", ratedpower=5E6)
        if method == "absolute"
            return first .- second
        elseif method == "normalizedindividually"
            absolute_errors = errors(first, second, method="absolute")
            normalized_errors = absolute_errors./first
            return normalized_errors
        elseif method == "normbyrated"
            absolute_errors = errors(first, second, method="absolute")
            normalized_errors = absolute_errors./ratedpower
            return normalized_errors
        elseif method == "normbyfirst"
            absolute_errors = errors(first, second, method="absolute")
            normalized_errors = absolute_errors./maximum(first)
            return normalized_errors
        elseif method == "normbyrow"
            absolute_errors = errors(first, second, method="absolute")
            normalized_errors = zeros(size(absolute_errors))
            for i in 1:length(first[:,1])
                normalized_errors[i,:] = absolute_errors[i,:]./maximum(absolute_errors[i,:])
            end
            return normalized_errors
        elseif method == "normbyselfrow"
            firstnormed = zeros(size(first))
            secondnormed = zeros(size(second))
            for i in 1:length(first[:,1])
                firstnormed[i,:] = firstnormed[i,:]/maximum(firstnormed[i,:])
                secondnormed[i,:] = secondnormed[i,:]/maximum(secondnormed[i,:])
            end
        
            normalized_errors = errors(firstnormed, secondnormed, method="absolute")
        
            return normalized_errors
        end
    end


    # load wake data 
    wakecountfile = "image-data/power/turbine-wakes-$case-$tuning.txt"
    
    # flowfarm base data 
    basepowerfileff = "image-data/power/turbine-power-ff-100pts-$case-$tuning.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/turbine-power-$case-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/turbine-power-ff-100pts-low-ti-$tuning-opt.txt"
    
    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/turbine-power-low-ti-les-opt.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    wake_count = transpose(readdlm(wakecountfile, ',', skipstart=1, header=false))
    baseff = transpose(readdlm(basepowerfileff, ',', skipstart=1, header=false))
        # DataFrame(CSV.File(basepowerfileff, datarow=2, header=false))
    optff = transpose(readdlm(optpowerfileff, ',', skipstart=1, header=false))
    basesowfa = transpose(readdlm(basepowerfilesowfa, skipstart=1, header=false))
    # DataFrame(CSV.File(basepowerfilesowfa, datarow=2, header=false))
    optsowfa = transpose(readdlm(optpowerfilesowfa, skipstart=1, header=false))
    winddf = DataFrame(CSV.File(winddatafile, datarow=2, header=false))

    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    # set vmin and vmax 
    vmax = 26
    vmin = -vmax

    # calculate turbine errors for base case
    turberror = errors(basesowfa, baseff, method="normbyfirst")
    ers = (sum(basesowfa[i,:] for i=1:12 ) .- sum(baseff[j,:] for j=1:12))
    # ers = (basesowfa .- baseff)
    println("base turbine error sum $case $tuning: $(sqrt(sum(ers.^2)))")
    # println(basesowfa)
    # println(baseff)
    data = convert.(Int64, round.(turberror.*100, digits=0))
    
    for i = 1:12
        sumterm = 0
        for j = 1:38
            if wake_count[i,j] > 0
                sumterm += data[i,j]
            end
        end
        # println(sumterm)
    end
    # plot error on heatmap
    plot_turbine_heatmap(data, winddf.d, vmin, vmax, wake_count=wake_count)

    if savefigs
        plt.savefig(image_directory*image_name*"-"*case*"-"*tuning*".pdf", transparent=true)
    end

    # calculate turbine errors for opt case
    turberror = errors(optsowfa, optff, method="normbyfirst")
    println("opt turbine error sum $case $tuning: $(sum(turberror))")
    data = convert.(Int64, round.(turberror.*100, digits=0))

    # plot error on heatmap
    plot_turbine_heatmap(data, winddf.d, vmin, vmax)

    if savefigs
        plt.savefig(image_directory*image_name*"-"*case*"-"*tuning*"-opt.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end

end

function horns_rev_rows_verification_figure(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="horns-rev-rows", nsamplepoints=100)

    # load FLOWFarm and Niayifar data
    datafile = "image-data/verification/horns-rev-rows-$nsamplepoints-sample-points.txt"
    data = readdlm(datafile, skipstart=1)
    rows = data[:, 1]
    normalized_power_les_niayifar = data[:,2]
    normalized_power_model_niayifar = data[:, 3]
    normalized_power_averaged_ff_no_ti = data[:,4]
    normalized_power_averaged_ff_ti = data[:, 5]

    fig, ax = plt.subplots(figsize=(5,3))

    ax.scatter(rows, normalized_power_les_niayifar, c=colors[1], label="Niayifar 2016 LES", marker="o")
    ax.scatter(rows, normalized_power_model_niayifar, c=colors[4], label="Niayifar 2016 Model", marker="*")
    ax.scatter(rows, normalized_power_averaged_ff_no_ti, edgecolors=colors[2], label="BP w/o Local TI", marker="^", facecolors="none")
    ax.scatter(rows, normalized_power_averaged_ff_ti, edgecolors=colors[3], label="BP w/Local TI", marker="v", facecolors="none")
    
    # format the figure
    ax.set(xlabel="Row", ylabel="Normalized Power", ylim=[0,1.1], xticks=Int.(rows))
    ax.legend(frameon=false,ncol=1)

    # remove upper and right bounding box
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    plt.tight_layout()

    if savefigs
        plt.savefig(image_directory*image_name*"-$nsamplepoints-sample-points.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function horns_rev_direction_verification_figure(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="horns-rev-directions", nsamplepoints=100)

    # load FLOWFarm and Niayifar data
    lesfile = "../src/inputfiles/results-niayifar-2016/horns-rev-normalized-power-by-direction-les.txt"
    modelfile = "image-data/verification/horns-rev-direction-$nsamplepoints-sample-points.txt"

    lesdata = readdlm(lesfile, ',', skipstart=1)
    modeldata = readdlm(modelfile, skipstart=1)
    
    lesdirections = lesdata[:, 1]
    normalized_power_les_niayifar = lesdata[:,2]

    modeldirections = modeldata[:, 1]
    normalized_power_model_niayifar = modeldata[:, 2]
    normalized_power_averaged_ff_no_ti = modeldata[:,3]
    normalized_power_averaged_ff_ti = modeldata[:, 4]

    fig, ax = plt.subplots(figsize=(5,3))
    println(lesdirections)
    println(modeldirections)
    ax.plot(lesdirections, normalized_power_les_niayifar, c=colors[1], label="Niayifar 2016 LES", marker="o", linestyle="none", markersize=4)
    ax.plot(modeldirections, normalized_power_model_niayifar, "-", c=colors[4], label="Niayifar 2016 Model")
    ax.plot(modeldirections, normalized_power_averaged_ff_no_ti, "--", c=colors[2], label="BP w/o Local TI")
    ax.plot(modeldirections, normalized_power_averaged_ff_ti, ":", c=colors[3], label="BP w/Local TI")
    
    # format the figure
    ax.set(xlabel="Direction (deg.)", ylabel="Normalized Power", ylim=[0,1])
    ax.legend(frameon=false,ncol=2)

    # remove upper and right bounding box
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    plt.tight_layout()

    if savefigs
        plt.savefig(image_directory*image_name*"-$nsamplepoints-sample-points.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end


function plot_circle(cx,cy,R,color,ax;fill=false,alpha=1.0,linestyle="-",linewidth=1,label="nolabel")
    N = 1000
    x = zeros(N)
    y = zeros(N)
    theta = range(0.0,stop=2*pi,length=N)
    for i=1:N
        x[i] = cx + cos(theta[i])*R
        y[i] = cy + sin(theta[i])*R
    end

    if fill==true
        if label=="nolabel"
            ax.fill(x,y,color=color,alpha=alpha, edgecolor="None")
        else
            ax.fill(x,y,color=color,alpha=alpha, edgecolor="None",label=label)
        end
    else
        if label=="nolabel"
            ax.plot(x,y,color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth)
        else
            ax.plot(x,y,color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth,label=label)
        end
    end
end

function rotor_sample_points(x1,y1,x2,y2;color="C0",fontsize=10,filename="nosave")

    plt.figure(figsize=(6,3))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    ax1.spines["right"].set_visible(false)
    ax1.spines["top"].set_visible(false)
    ax2.spines["right"].set_visible(false)
    ax2.spines["top"].set_visible(false)

    # ax1.plot(x1,y1,"o",)
    div = 130.0
    bladeX = [3.,7.,10.,15.,20.,25.,30.,35.,30.,25.,20.,15.,10.,5.,3.,3.]
    bladeY = [0.,0.,0.8,1.5,1.7,1.9,2.1,2.3,2.4,2.4,2.4,2.4,2.4,2.4,2.4,0.].-1.5
    #
    T = 110/div
    H = 0/div
    D = 130/div
    R = 65/div
    d = [6.3,5.5,4.]./div

    c1 = R/35

    px1 = [0.0-d[1]/2,0.0-d[2]/2,0.0-d[3]/2,0.0+d[3]/2,0.0+d[2]/2,0.0+d[1]/2,0.0-d[1]/2]
    py1 = [-T,-T/2,-3.0*c1,-3.0*c1,-T/2,-T,-T]
    #
    #add blades

    #
    angle1 = 5.

    blade1X = bladeX*cos(deg2rad(angle1))-bladeY*sin(deg2rad(angle1))
    blade1Y = bladeX*sin(deg2rad(angle1))+bladeY*cos(deg2rad(angle1))

    blade2X = bladeX*cos(deg2rad(angle1+120.))-bladeY*sin(deg2rad(angle1+120.))
    blade2Y = bladeX*sin(deg2rad(angle1+120.))+bladeY*cos(deg2rad(angle1+120.))

    blade3X = bladeX*cos(deg2rad(angle1+240.))-bladeY*sin(deg2rad(angle1+240.))
    blade3Y = bladeX*sin(deg2rad(angle1+240.))+bladeY*cos(deg2rad(angle1+240.))
    #
    plot_circle(0,H,R,color,ax1,alpha=0.5)
    ax1.plot(px1,py1,color=color, linewidth=1.2*1.5,alpha=0.5)
    plot_circle(0,H,3*c1,color,ax1,alpha=0.5)
    ax1.plot(blade1X.*c1.+0.0, blade1Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax1.plot(blade2X.*c1.+0.0, blade2Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax1.plot(blade3X.*c1.+0.0, blade3Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax1.axis("equal")

    angle2 = -20.

    blade1X = bladeX*cos(deg2rad(angle2))-bladeY*sin(deg2rad(angle2))
    blade1Y = bladeX*sin(deg2rad(angle2))+bladeY*cos(deg2rad(angle2))

    blade2X = bladeX*cos(deg2rad(angle2+120.))-bladeY*sin(deg2rad(angle2+120.))
    blade2Y = bladeX*sin(deg2rad(angle2+120.))+bladeY*cos(deg2rad(angle2+120.))

    blade3X = bladeX*cos(deg2rad(angle2+240.))-bladeY*sin(deg2rad(angle2+240.))
    blade3Y = bladeX*sin(deg2rad(angle2+240.))+bladeY*cos(deg2rad(angle2+240.))
    #
    plot_circle(0,H,R,color,ax2,alpha=0.5)
    ax2.plot(px1,py1,color=color, linewidth=1.2*1.5,alpha=0.5)
    plot_circle(0,H,3*c1,color,ax2,alpha=0.5)
    ax2.plot(blade1X.*c1.+0.0, blade1Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax2.plot(blade2X.*c1.+0.0, blade2Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax2.plot(blade3X.*c1.+0.0, blade3Y.*c1.+H, linewidth=1*1.5, color=color,alpha=0.5)
    ax2.axis("equal")


    ax1.plot(x1,y1,"o",color=color,markersize=4,label="sampling points")
    ax2.plot(x2,y2,"o",color=color,markersize=4)
    ax1.legend(fontsize=fontsize)

    # ax1.set_xlabel(r"Horizontal Distance From Hub, $\Delta y/d$",fontsize=fontsize)
    ax1.set_xlabel("Horizontal Distance From Hub, "*L"\Delta y/d",fontsize=fontsize)
    ax1.set_ylabel("Vertical Distance From Hub, "*L"\Delta z/d",fontsize=fontsize)
    ax2.set_xlabel("Horizontal Distance From Hub, "*L"\Delta y/d",fontsize=fontsize)
    ax2.set_ylabel("Vertical Distance From Hub, "*L"\Delta z/d",fontsize=fontsize)

    ax1.text(0,-1.3,"(a)",horizontalalignment="center",fontsize=fontsize)
    ax2.text(0,-1.3,"(b)",horizontalalignment="center",fontsize=fontsize)

    plt.subplots_adjust(top=0.99,bottom=0.22,right=0.92,left=0.12,wspace=0.8)

    if filename != "nosave"
        plt.savefig(filename,transparent=true)
    end
end

function sunflower_points(n, alpha=1.0)
    # this function generates n points within a circle in a sunflower seed pattern
    # the code is based on the example found at
    # https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

    function radius(k, n, b)
        if (k + 1) > n - b
            r = 1.0 # put on the boundary
        else
            r = sqrt((k + 1.0) - 1.0 / 2.0) / sqrt(n - (b + 1.0) / 2.0)  # apply squareroot
        end
        return r
    end


    x = zeros(n)
    y = zeros(n)

    b = round(alpha * sqrt(n)) # number of boundary points

    phi = (sqrt(5.0) + 1.0) / 2.0  # golden ratio

    for k = 1:n
        r = radius(k, n, b)

        theta = 2.0 * pi * (k+1) / phi^2

        x[k] = r * cos(theta)
        y[k] = r * sin(theta)
    end

    return x, y
end

function turbine_layouts(colors ;rotor_diameter=126.4,les_side=5000,
                                    c1="C0",c2="C1",color="C0",fontsize=10,
                                    case="low-ti-opt", showfigs=false, savefigs=false)

    les_radius = les_side/2.0
    boundary_radius = les_radius - 500.0 - rotor_diameter/2.0

    # load data
    if case == "low-ti-opt"
        datafile = "image-data/layouts/opt/optresultsmilestone.csv"
        data = readdlm(datafile, ',', skipstart=1).+les_radius
    elseif case == "low-ti"
        datafile = "../src/inputfiles/farms/layout_38turb_round.txt"
        data = readdlm(datafile, skipstart=1).*rotor_diameter .+ (les_radius - boundary_radius)
    else
        ErrorException("Case not available")
    end

    xlocs = data[:,1]
    ylocs = data[:,2]

    plt.figure(figsize=(5,3))
    ax = plt.gca()

    ax.axes.xaxis.set_visible(false)
    ax.axes.yaxis.set_visible(false)

    ax.spines["top"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.spines["left"].set_visible(false)
    ax.spines["right"].set_visible(false)

    # xaxis.set_visible(False)

    R = rotor_diameter/2
    cx = 2500
    cy = 2500
    lw = 0.75
    plot_circle(cx,cy,boundary_radius,colors[2],ax,linestyle="--",linewidth=lw,label="Farm boundary")
    ax.annotate("Farm Boundary", (les_radius, les_radius-boundary_radius-rotor_diameter*2.5), color=colors[2])

    les_x = [0,les_side,les_side,0,0]
    les_y = [0,0,les_side,les_side,0]
    ax.plot(les_x,les_y,"-",color=colors[4],linewidth=lw,label="LES Domain")
    ax.annotate("LES Domain", (0, les_side+rotor_diameter/2.0), color=colors[4])
    
    nturbs = length(xlocs)
    for i=1:nturbs
        plot_circle(xlocs[i],ylocs[i],R,colors[3],ax,linewidth=1,fill=false)
        ax.text(xlocs[i]+rotor_diameter/5,ylocs[i]+rotor_diameter/5,"$i",
                fontsize=fontsize-2,horizontalalignment="left",verticalalignment="bottom",
                color=colors[1])
    end

    # ax.legend(fontsize=fontsize,bbox_to_anchor=[1.0, 1.0])
    ax.axis("square")
    # ax.set_xlim(-100,les_side+1000)

    # ax.set_xlabel("Turbine X Position, m",fontsize=fontsize)
    # ax.set_ylabel("Turbine Y Position, m",fontsize=fontsize)
    
    plt.subplots_adjust(left=0.15,bottom=0.12,top=0.99,right=0.6)

    plt.tight_layout()

    if savefigs
        plt.savefig("images/"*case*"-layout.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function windrose(d1,f1,d2,f2;color="C0",alpha=0.5,fontsize=8,filename="nosave")
    
    f1 = f1./sum(f1)
    f2 = f2./sum(f2)
    if maximum(d1) > 100
        d1 = deg2rad.(d1)
    end
    if maximum(d2) > 100
        d2 = deg2rad.(d2)
    end
    println(wd1)
    plt.figure(figsize=(6,3))
    ax1 = plt.subplot(121,projection="polar")
    ax2 = plt.subplot(122,projection="polar")

    ndirs1 = length(d1)
    width1 = 2*pi/ndirs1
    ax1.bar(pi/2 .-d1,f1,width=width1,color=color,alpha=alpha,edgecolor="black")
    ax1.set_xticks((0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
    ax1.set_xticklabels(("E","NE","N","NW","W","SW","S","SW"),fontsize=fontsize)
    ax1.set_rgrids((0.04,0.08,0.12),("4%","8%","12%"),angle=-20,fontsize=fontsize)
    for tick in ax1.yaxis.get_majorticklabels()
        tick.set_horizontalalignment("center")
    end


    ndirs2 = length(d2)
    width2 = 2*pi/ndirs2
    ax2.bar(pi/2 .-d2,f2,width=width2,color=color,alpha=alpha,edgecolor="black")
    ax2.set_xticks((0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
    ax2.set_xticklabels(("E","NE","N","NW","W","SW","S","SW"),fontsize=fontsize)
    ax2.set_rgrids((0.02,0.035,0.05),("2%","3.5%","5%"),angle=-20,fontsize=fontsize)
    for tick in ax2.yaxis.get_majorticklabels()
        tick.set_horizontalalignment("center")
    end

    ax1.set_title("(a)", y=-0.25,fontsize=fontsize)
    ax2.set_title("(b)", y=-0.25,fontsize=fontsize)

    plt.subplots_adjust(left=0.05,right=0.95,top=0.9,bottom=0.2)

    if filename != "nosave"
        plt.savefig(filename,transparent=true)
    end

end

function vertical_slice(colors; savefigs=false, showfigs=false, version="interpolated")

    # set input values 
    diam = 126.4

    # load flowfield velocities
    ffvelocities = readdlm("image-data/verification/vertical-slice-interpolated.txt")

    # load plot ranges
    ranges = readdlm("image-data/verification/vertical-slice-ranges.txt")
    xres = ranges[1, 2]
    zres = ranges[1, 3]
    xmax = ranges[2, 2]
    zmax = ranges[2, 3]
    xmin = ranges[3, 2]
    zmin = ranges[3, 3]

    # set up point grid for flow field
    xrange = xmin:(xmax-xmin)/xres:xmax
    zrange = zmin:(zmax-zmin)/zres:zmax

    # create meshgrid 
    xg, zg = ml.meshgrid(collect(xrange), collect(zrange))

    # set contour levels
    # levels = 1:8.2
    levels = 0:0.5:9

    # get colormap 
    cmap = custum_color_map(idx=[3,2])
    
    # generate contour plot
    fig, ax = plt.subplots(figsize=(7,2))
    
    
    cs = ax.contourf(xg./diam, zg./diam, ffvelocities, levels, cmap=cmap)

    if version == "original"
        r = plt.matplotlib.patches.Rectangle((0,0),725.4980208718556/diam,2, color="w")
        ax.add_patch(r)
    end

    cbar = ax.figure.colorbar(cs, ax=ax, label="Wind Speed (m/s)", orientation="vertical")
    # cbar.ax.set_xlabel("Wind Speed (m/s)", rotation=270)
    
    println(minimum(ffvelocities))
    println(cs.levels)
    # ax.clabel(cs, 4:8, inline=1) 

    # add turbine 
    radius = diam/2
    ax.plot([0.0,0.0], [90-radius, 90+radius]./diam, linewidth=5, color="k")

    # format figure 
    ax.set(xticks=0:4:20, yticks=0:2)
    ax.set(xlabel=L"$x/d_0$", ylabel=L"$z/d_0$")

    plt.tight_layout()

    if savefigs
        plt.savefig("images/vertical-slice-$version.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function generate_images_for_publication()

    fontsize = 8
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]

    rcParams = PyPlot.matplotlib.rcParams
    rcParams["font.size"] = fontsize
    rcParams["lines.markersize"] = 8
    rcParams["axes.prop_cycle"] = colors

    savefigs = true 
    showfigs = true

    # wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti")
    # wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti")
    # layout(colors, fontsize)
    # opt_comparison_table()
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs)
    turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="alldirections")
    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="alldirections")
    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="all")
    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="all")
    # horns_rev_rows_verification_figure(colors, fontsize, nsamplepoints=1, savefigs=savefigs, showfigs=showfigs)
    # horns_rev_rows_verification_figure(colors, fontsize, nsamplepoints=100, savefigs=savefigs, showfigs=showfigs)
    # horns_rev_direction_verification_figure(colors, fontsize, nsamplepoints=1, savefigs=savefigs, showfigs=showfigs)
    # horns_rev_direction_verification_figure(colors, fontsize, nsamplepoints=100, savefigs=savefigs, showfigs=showfigs)

    # turbine_layouts(colors, case="low-ti-opt", showfigs=showfigs, savefigs=savefigs)
    # turbine_layouts(colors, case="low-ti", showfigs=showfigs, savefigs=savefigs)

    # vertical_slice(colors, savefigs=savefigs, showfigs=showfigs)
    # vertical_slice(colors, savefigs=savefigs, showfigs=showfigs, version = "original")
end