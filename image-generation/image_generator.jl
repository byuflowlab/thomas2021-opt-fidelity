import PyPlot; const plt = PyPlot
using DataFrames 
using CSV
using Colors, ColorSchemes
using FLOWFarm; const ff = FLOWFarm
using PrettyTables
using DelimitedFiles

function custum_color_map()
    colors = [colorant"#BDB8AD", colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    # @pyimport matplotlib.colors as matcolors
    # cmap = matcolors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],"A")

    return plt.ColorMap("BlueGrayOrange", [colors[3],colors[1],colors[4]])
end

function heatmap(data, row_labels, col_labels; ax=nothing, cbar_kw=Dict(), cbarlabel="", use_cbar=true, labelpixels=true, vcolor="w", edgecolor="w", boxmaxmin=true)
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
    rotor_diameter = 126.4 # m 
    hub_height = 90.0 # m 
    swept_top = hub_height + rotor_diameter/2.0
    swept_bottom = hub_height - rotor_diameter/2.0

    # create plot 
    fig, ax = plt.subplots(figsize=(10.5,4.5))

    # plot les data 
    ax.scatter(df_les.s, df_les.h, color=colors[2])

    if case=="high-ti"
        ax.annotate("LES", (7.75, 125), color=colors[2], alpha=1.0, size=fontsize)
    elseif case=="low-ti"
        ax.annotate("LES", (7.9, 125), color=colors[2], alpha=1.0, size=fontsize)
    end

    # plot model 
    ax.plot(df_model.s, df_model.h, color=colors[3])

    if case=="high-ti"
        ax.annotate("Curve Fit", (8.15, 95), color=colors[3], alpha=1.0, size=fontsize)
    elseif case=="low-ti"
        ax.annotate("Curve Fit", (8.25, 95), color=colors[3], alpha=1.0, size=fontsize) 
    end
    # add rotor swept area 
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_top, swept_top], color=colors[1], linestyle="--")
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_bottom, swept_bottom], color=colors[1], linestyle="--")
    ax.fill_between(5:9, swept_bottom,swept_top, alpha=0.15, color=colors[1])
    ax.annotate("Rotor-Swept Region", (5.5, (hub_height+swept_bottom)/2.0), color=colors[1], alpha=1.0, size=fontsize)

    # add hub height
    plt.plot([5,9], [hub_height, hub_height], color=colors[1], linestyle="-.")

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
    ax.set_xlim([5.0, 9.0])
    ax.set_ylim([0.0, 310.0])

    # set ticks 
    plt.xticks(5:9, size=fontsize)
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

function turbine_comparison_figures(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="turbine-comparison", case="low-ti")
    
    function plot_turbine_heatmap(data, winddirections, vmin, vmax)

        # number of turbines in the farm
        nturbines = 38

        # intialize figure
        fig, ax = plt.subplots(figsize=(8,4))

        # set tick locations and labels
        ticks = vmin:4:vmax

        println(vmax)
        
        # generate a custom color map 
        cmap = custum_color_map()

        # generate dict of colormap options
        d = Dict(:shrink => 0.47, :ticks=>ticks, :aspect=>20, :orientation=>"horizontal", :cmap=>cmap, :pad=>0.05)
        
        # set row labels
        rowlabels = convert.(Int64, round.((winddirections), digits=0))

        # create heatmap
        im, cbar = heatmap(data, rowlabels, 1:nturbines, ax=ax,
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
    vmax = 22
    vmin = -vmax

    # calculate turbine errors for base case
    turberror = errors(basesowfa, baseff, method="normbyfirst")
    data = convert.(Int64, round.(turberror.*100, digits=0))

    # plot error on heatmap
    plot_turbine_heatmap(data, winddf.d, vmin, vmax)

    if savefigs
        plt.savefig(image_directory*image_name*"-"*case*".pdf", transparent=true)
    end

    # calculate turbine errors for opt case
    turberror = errors(optsowfa, optff, method="normbyfirst")
    data = convert.(Int64, round.(turberror.*100, digits=0))

    # plot error on heatmap
    plot_turbine_heatmap(data, winddf.d, vmin, vmax)

    if savefigs
        plt.savefig(image_directory*image_name*"-"*case*"-opt.pdf", transparent=true)
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
    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs)
    horns_rev_rows_verification_figure(colors, fontsize, nsamplepoints=1, savefigs=savefigs, showfigs=showfigs)
    horns_rev_rows_verification_figure(colors, fontsize, nsamplepoints=100, savefigs=savefigs, showfigs=showfigs)
    horns_rev_direction_verification_figure(colors, fontsize, nsamplepoints=1, savefigs=savefigs, showfigs=showfigs)
    horns_rev_direction_verification_figure(colors, fontsize, nsamplepoints=100, savefigs=savefigs, showfigs=showfigs)
end