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

function wind_shear_tuning(colors, fontsize; ax=nothing, showfigs=false, savefigs=false, image_directory="images/", image_name="windshear-", case="high-ti")

    # load data 
    df_model = DataFrame(CSV.File("./image-data/wind-shear/$(case)/wind-shear-tuned-$(case).txt", skipto=2, header=false))
    df_les = DataFrame(CSV.File("./image-data/wind-shear/$(case)/wind-shear-les-$(case).txt", skipto=2, header=false))

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
    if ax === nothing
        fig, ax = plt.subplots(figsize=(5.5,2.5))
    end
    # plot les data 
    ax.scatter(df_les.s, df_les.h, color=colors[2], s=10)

    case=="high-ti" && ax.annotate("High-TI LES", (8.8, 130), color=colors[2], alpha=1.0, size=fontsize)
    case=="low-ti" && ax.annotate("Low-TI LES", (8.425, 130), color=colors[2], alpha=1.0, size=fontsize)

    # plot model 
    ax.plot(df_model.s, df_model.h, color=colors[3])

    case=="high-ti" && ax.annotate("High-TI Curve Fit", (6.7, 130), color=colors[3], alpha=1.0, size=fontsize)
    case=="low-ti" && ax.annotate("Low-TI Curve Fit", (6.5, 130), color=colors[3], alpha=1.0, size=fontsize)
    
    # add rotor swept area 
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_top, swept_top], color=colors[1], linestyle="--")
    # ax.plot([minimum(df_model.s), maximum(df_model.s)], [swept_bottom, swept_bottom], color=colors[1], linestyle="--")
    ax.fill_between(minv:maxv, swept_bottom,swept_top, alpha=0.15, color=colors[1])
    ax.annotate("Rotor-Swept Region", (4.5, (hub_height+swept_bottom)/2.0), color=colors[1], alpha=1.0, size=fontsize)

    # add hub height
    ax.plot([minv,maxv], [hub_height, hub_height], color=colors[1], linestyle="-.")

    # add axis labels
    ax.set_xlabel(L"Wind speed (m s$^{-1}$)")#, fontsize=fontsize)
    ax.set_ylabel("Height (m)")#, fontsize=fontsize)

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
    ax.set(xticks=minv:maxv)
    ax.set(yticks=0:100:300)

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

function wind_shear_tuning_compound(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="windshear-compound")
    fig, ax = plt.subplots(1,2, figsize=(10,3))

    wind_shear_tuning(colors, fontsize*1.25, ax=ax[1], case="high-ti")
    wind_shear_tuning(colors, fontsize*1.25, ax=ax[2], case="low-ti")

    ax[1].text(7,-120,"(a)",horizontalalignment="center",fontsize=fontsize*1.25)
    ax[2].text(7,-120,"(b)",horizontalalignment="center",fontsize=fontsize*1.25)
    
    plt.tight_layout()

    # save figure
    if showfigs
        plt.show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*".pdf", transparent=true)
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
    df = DataFrame(CSV.File(datafile, skipto=2, header=false))

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

function opt_comparison_table(;case="high-ti", tuning="sowfa-nrel")

    if case == "high-ti"
        layoutid = 83
        n = 4
    elseif case == "low-ti"
        layoutid = 252
        n = 2
    end

    # flowfarm base data 
    basepowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout1-base.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout$layoutid-opt$n.txt"

    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les-opt$n.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    baseffdf = DataFrame(CSV.File(basepowerfileff, skipto=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, skipto=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, skipto=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, skipto=2, header=false))
    winddf = DataFrame(CSV.File(winddatafile, skipto=2, header=false))

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
        println("$bs MW &")
        println("$bf MW &")
        println("$bc \$\\%\$ & &")
        println("$os MW &")
        println("$of MW &")
        println("$oc \$\\%\$ & &")
        println("$is \$\\%\$ &")
        println("$iff \$\\%\$ \\\\")
        
        # println("\\multicolumn{1}{l}{\$$(Int(round(winddf.d[i], digits=0)))^{\\circ}\$} & &")
        # println("\\SI[per-mode=symbol]{$bs}{\\mega\\watt} &")
        # println("\\SI[per-mode=symbol]{$bf}{\\mega\\watt} &")
        # println("\\SI[per-mode=symbol]{$bc}{\\percent} & &")
        # println("\\SI[per-mode=symbol]{$os}{\\mega\\watt} &")
        # println("\\SI[per-mode=symbol]{$of}{\\mega\\watt} &")
        # println("\\SI[per-mode=symbol]{$oc}{\\percent} & &")
        # println("\\SI[per-mode=symbol]{$is}{\\percent} &")
        # println("\\SI[per-mode=symbol]{$iff}{\\percent} \\\\")
        
        # println(dirdata[i,:])
    end

    println("\\multicolumn{1}{l}{AEP} & &")
    println("$(round(aepbasesowfa*1e-9, digits=1)) GW h &")
    println("$(round(aepbaseff*1e-9, digits=1)) GW h &")
    println("$(round(aepdiffbase, digits=1)) \$\\%\$ & &")
    println("$(round(aepoptsowfa*1e-9, digits=1)) GW h &")
    println("$(round(aepoptff*1e-9, digits=1)) GW h &")
    println("$(round(aepdiffopt, digits=1)) \$\\%\$ & &")
    println("$(round(aepimpsowfa, digits=1)) \$\\%\$ & ")
    println("$(round(aepimpff, digits=1)) \$\\%\$ \\\\")
    
    # println("\\multicolumn{1}{l}{AEP} & &")
    # println("\\SI[per-mode=symbol]{$(round(aepbasesowfa*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    # println("\\SI[per-mode=symbol]{$(round(aepbaseff*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    # println("\\SI[per-mode=symbol]{$(round(aepdiffbase, digits=1))}{\\percent} & &")
    # println("\\SI[per-mode=symbol]{$(round(aepoptsowfa*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    # println("\\SI[per-mode=symbol]{$(round(aepoptff*1e-9, digits=1))}{\\giga\\watt\\hour} &")
    # println("\\SI[per-mode=symbol]{$(round(aepdiffopt, digits=1))}{\\percent} & &")
    # println("\\SI[per-mode=symbol]{$(round(aepimpsowfa, digits=1))}{\\percent} &")
    # println("\\SI[per-mode=symbol]{$(round(aepimpff, digits=1))}{\\percent} \\\\")

end

function directional_comparison_figure(colors, fontsize; ax = nothing, showfigs=false, savefigs=false, image_directory="images/", image_name="directional-comparison", case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="power")

    if case == "low-ti"
        layout = 252
        n = 2
        # layout = 385
        # n = 3
        sn = 2
        
    elseif case == "high-ti"
        layout = 83
        n = 4
        sn = 4
    end
    image_name="directional-comparison-$case-$plottype"
    # load wake data 
    wakecountfilebase = "image-data/wakes/turbine-wakes-base-ff-$case-$tuning-layout1.txt"
    wakecountfileopt = "image-data/wakes/turbine-wakes-opt$n-ff-$case-$tuning-layout$layout.txt"

    # flowfarm base data 
    basepowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout1-base.txt"
    # basepowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout1-base-8ms.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout$layout-opt$n.txt"

    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les-opt$sn.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    basewcdf = DataFrame(CSV.File(wakecountfilebase, skipto=1, header=false))
    optwcdf = DataFrame(CSV.File(wakecountfileopt, skipto=1, header=false))
    baseffdf = DataFrame(CSV.File(basepowerfileff, skipto=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, skipto=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, skipto=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, skipto=2, header=false))
    winddf = DataFrame(CSV.File(winddatafile, skipto=2, header=false))

    println("compare")
    println(sum.(eachcol(baseffdf .- basesowfadf)))
    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    nturbines = length(basesowfadf.Column1)
    println("nturbs = ", nturbines)

    ndirections = length(winddf.d)

    # compute directional data 
    if normalize
        normalization_factor_sowfa = nturbines*maximum(maximum(eachrow(basesowfadf)))
        normalization_factor_bp = nturbines*maximum(maximum(eachrow(baseffdf)))
    else
        normalization_factor_sowfa = 1E6
        normalization_factor_bp = 1E6
    end

    println(normalization_factor_sowfa)
    println(normalization_factor_bp)

    basedirpowerff = sum.(eachcol(baseffdf))./normalization_factor_bp
    optdirpowerff = sum.(eachcol(optffdf))./normalization_factor_bp
    basedirpowersowfa = sum.(eachcol(basesowfadf))./normalization_factor_sowfa
    optdirpowersowfa = sum.(eachcol(optsowfadf))./normalization_factor_sowfa

    # create directional power bar charts
    compound = true
    if ax === nothing
        fig, ax = plt.subplots(1, figsize=[6,4])
        compound=false
    end
    
    if plottype == "power"
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , optdirpowersowfa, color=colors[3], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("SOWFA-optimized", (160, 65), color=colors[3]))
        compound || (case=="low-ti" && ax.annotate("SOWFA-optimized", (160, 62), color=colors[3]))

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , basedirpowersowfa, color=colors[2], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("SOWFA-base", (160, 55.5), color=colors[2]))
        compound || (case=="low-ti" && ax.annotate("SOWFA-base", (160, 54.5), color=colors[2]))
        
        # ax.bar(winddf.d .+ 5, optdirpowerff.*1E-6, label="BP-Opt", width=10, color=colors[4])
        ax.plot(winddf.d , optdirpowerff, color=colors[4], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("BP-optimized", (160, 60), color=colors[4]))
        compound || (case=="low-ti" && ax.annotate("BP-optimized", (160, 56.5), color=colors[4]))

        # ax.bar(winddf.d .+ 5, basedirpowerff.*1E-6, label="BP-Base", width=10, color=colors[1])
        ax.plot(winddf.d , basedirpowerff, color=colors[1], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("BP-base", (160, 53), color=colors[1]))
        compound || (case=="low-ti" && ax.annotate("BP-base", (160, 50), color=colors[1]))

        # format the figure
        compound || ax.set(xticks=winddf.d, ylim=[50, 70], xlabel="Direction (degrees)", ylabel="Directional power (MW)")
        # ax.legend(frameon=false,ncol=2)
    elseif plottype == "improvement"
        
        improvementsowfa = (optdirpowersowfa./basedirpowersowfa .- 1.0)*100
        improvementff = (optdirpowerff./basedirpowerff .- 1.0)*100
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , improvementsowfa, color=colors[3], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("SOWFA", (100, 9), color=colors[3]))
        compound || (case=="low-ti" && ax.annotate("SOWFA", (100, 10), color=colors[3]))

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , improvementff, color=colors[2], marker="o", linestyle="--")
        compound || (case=="high-ti" && ax.annotate("BP", (160, 6), color=colors[2]))
        compound || (case=="low-ti" && ax.annotate("BP", (160, 9), color=colors[2]))

        # format the figure
        compound || ax.set(xticks=winddf.d, ylim=[0, 15], xlabel="Direction (degrees)", ylabel="Directional Improvment (%)")

    elseif plottype == "error"
        errorbase = (basedirpowerff./basedirpowersowfa .- 1.0)*100
        erroropt = (optdirpowerff./optdirpowersowfa .- 1.0)*100
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , errorbase, color=colors[3], marker="o", linestyle="--", label="Base")
        compound || (case=="high-ti" && ax.annotate("Base", (130, -4.5), color=colors[3]))
        compound || (case=="low-ti" && ax.annotate("Base", (130, -0.5), color=colors[3]))

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , erroropt, color=colors[2], marker="o", linestyle="--", label="Optimized")
        compound || (case=="high-ti" && ax.annotate("Optimized", (130, -8), color=colors[2]))
        compound || (case=="low-ti" && ax.annotate("Optimized", (130, -2.6), color=colors[2]))

        # format the figure
        compound || ax.set(xticks=winddf.d, ylim=[-10, 2], xlabel="Direction (degrees)", ylabel="Directional Error (%)")
        
    elseif plottype == "wakeloss" || plottype == "annualenergyloss" || plottype == "annualenergylossincrease"

        # find which turbines are not waked
        # baseffarray = convert(Array{Float64, (nturbines, ndirections)}, baseffdf)
        # optwcarray = convert(Array, optwcdf)
        ideal_power_base_ff = zeros(ndirections)
        ideal_power_opt_ff = zeros(ndirections)
        ideal_power_base_sowfa = zeros(ndirections)
        ideal_power_opt_sowfa = zeros(ndirections)

        # find ideal power
        for j = eachindex(winddf.d)
            nb = 0
            no = 0
            for i = 1:nturbines
                if basewcdf[i,j] == 0.0
                    ideal_power_base_ff[j] += baseffdf[i,j]./normalization_factor_bp 
                    ideal_power_base_sowfa[j] += basesowfadf[i,j]./normalization_factor_sowfa
                    nb += 1
                end
                if optwcdf[i,j] == 0.0
                    ideal_power_opt_ff[j] += optffdf[i,j]./normalization_factor_bp
                    ideal_power_opt_sowfa[j] += optsowfadf[i,j]./normalization_factor_sowfa
                    no += 1
                end
            end
            ideal_power_base_ff[j] /= nb
            ideal_power_base_sowfa[j] /= nb
            ideal_power_opt_ff[j] /= no
            ideal_power_opt_sowfa[j] /= no
        end
        
        ideal_power_base_ff .*= nturbines
        ideal_power_base_sowfa .*= nturbines
        ideal_power_opt_ff .*= nturbines
        ideal_power_opt_sowfa .*= nturbines

        # get wake loss 
        basewakelossff = 100.0.*(1.0 .- basedirpowerff./ideal_power_base_ff)
        basewakelosssowfa = 100.0.*(1.0 .- basedirpowersowfa./ideal_power_base_sowfa)
        optwakelossff = 100.0.*(1.0 .- optdirpowerff./ideal_power_opt_ff)
        optwakelosssowfa = 100.0.*(1.0 .- optdirpowersowfa./ideal_power_opt_sowfa)

        if plottype == "annualenergyloss"
            ts = winddf.p*365*24.0.*1E-3.*(1/100)
            basewakelossff .*= ts.*ideal_power_base_ff
            basewakelosssowfa .*= ts.*ideal_power_base_sowfa
            optwakelossff .*= ts.*ideal_power_opt_ff
            optwakelosssowfa .*= ts.*ideal_power_opt_sowfa
            ts = ts[6]
        elseif plottype == "wakeloss"
            ts = 1
        elseif plottype == "annualenergylossincrease"
            ts = winddf.p*365*24.0.*1E-3.*(1/100)
            basewakelossff .*= ts.*ideal_power_base_ff
            basewakelosssowfa .*= ts.*ideal_power_base_sowfa
            optwakelossff .*= ts.*ideal_power_opt_ff
            optwakelosssowfa .*= ts.*ideal_power_opt_sowfa
            increaseff = -(optwakelossff - basewakelossff)
            increasesowfa = -(optwakelosssowfa - basewakelosssowfa)
            ts = ts[6]
        end

        # plot
        if plottype == "annualenergylossincrease"
            ax.plot(winddf.d , increasesowfa, color=colors[3], marker="o", linestyle="--", label="SOWFA")
            ax.plot(winddf.d , increaseff, color=colors[2], marker="o", linestyle="--", label="BP")
        else
            ax.plot(winddf.d , optwakelosssowfa, color=colors[3], marker="o", linestyle="--", label="SOWFA-optimized")
            if plottype == "annualenergyloss"
                compound || (case=="high-ti" && ax.annotate("SOWFA-optimized", (210, 4.5), color=colors[3]))
                compound || (case=="low-ti" && ax.annotate("SOWFA-optimized", (200, 6), color=colors[3]))
            elseif plottype == "wakeloss"
                compound || (case=="high-ti" && ax.annotate("SOWFA-optimized", (160, 5*ts), color=colors[3]))
                compound || (case=="low-ti" && ax.annotate("SOWFA-optimized", (160, 10*ts), color=colors[3]))
            elseif plottype == "annualenergylossincrease"

            end
            # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
            ax.plot(winddf.d , basewakelosssowfa, color=colors[2], marker="o", linestyle="--", label="SOWFA-base")
            if plottype == "annualenergyloss"
                compound || (case=="high-ti" && ax.annotate("SOWFA-base", (196, 10), color=colors[2]))
                compound || (case=="low-ti" && ax.annotate("SOWFA-base", (280, 9.5), color=colors[2]))
            else
                compound || (case=="high-ti" && ax.annotate("SOWFA-base", (160, 14.5*ts), color=colors[2]))
                compound || (case=="low-ti" && ax.annotate("SOWFA-base", (160, 19*ts), color=colors[2]))
            end

            # ax.bar(winddf.d .+ 5, optdirpowerff.*1E-6, label="BP-Opt", width=10, color=colors[4])
            ax.plot(winddf.d , optwakelossff, color=colors[4], marker="o", linestyle="--", label="BP-optimized")
            if plottype == "annualenergyloss"
                compound || (case=="high-ti" && ax.annotate("BP-optimized", (190, 8.5), color=colors[4]))
                compound || (case=="low-ti" && ax.annotate("BP-optimized", (195, 13), color=colors[4]))
            else
                compound || (case=="high-ti" && ax.annotate("BP-optimized", (160, 12.75*ts), color=colors[4]))
                compound || (case=="low-ti" && ax.annotate("BP-optimized", (160, 17.0*ts), color=colors[4]))
            end

            # ax.bar(winddf.d .+ 5, basedirpowerff.*1E-6, label="BP-Base", width=10, color=colors[1])
            ax.plot(winddf.d , basewakelossff, color=colors[1], marker="o", linestyle="--", label="BP-base")
            if plottype == "annualenergyloss"
                compound || (case=="high-ti" && ax.annotate("BP-base", (256, 12.3), color=colors[1]))
                compound || (case=="low-ti" && ax.annotate("BP-base", (265, 15), color=colors[1]))
            else
                compound || (case=="high-ti" && ax.annotate("BP-base", (160, 20*ts), color=colors[1]))
                compound || (case=="low-ti" && ax.annotate("BP-base", (160, 25.5*ts), color=colors[1]))
            end

            # format the figure
            if plottype == "annualenergyloss"
                compound || ax.set(xticks=winddf.d[1:2:end], yticks=collect(0:5:20), ylim=[0, 20], ylabel="Annual Directional Energy Loss (GW h)")
            else
                compound || ax.set(xticks=winddf.d[1:2:end], yticks=collect(0:5:30), ylim=[0, 30], ylabel="Directional Wake Loss (%)")
            end
            if !compound 
                ax.set(xlabel="Direction (degrees)")
            end
        end
    end

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

function directional_comparison_figure_polar(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="directional-comparison", case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="power")

    if case == "low-ti"
        layout = 252
        n = 2
    elseif case == "high-ti"
        layout = 83
        n = 4
    end

    # load wake data 
    wakecountfilebase = "image-data/wakes/turbine-wakes-base-ff-$case-$tuning-layout1.txt"
    wakecountfileopt = "image-data/wakes/turbine-wakes-opt$n-ff-$case-$tuning-layout$layout.txt"

    # flowfarm base data 
    basepowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout1-base.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout$layout-opt$n.txt"

    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les-opt$n.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    basewcdf = DataFrame(CSV.File(wakecountfilebase, skipto=1, header=false))
    optwcdf = DataFrame(CSV.File(wakecountfileopt, skipto=1, header=false))
    baseffdf = DataFrame(CSV.File(basepowerfileff, skipto=2, header=false))
    optffdf = DataFrame(CSV.File(optpowerfileff, skipto=2, header=false))
    basesowfadf = DataFrame(CSV.File(basepowerfilesowfa, skipto=2, header=false))
    optsowfadf = DataFrame(CSV.File(optpowerfilesowfa, skipto=2, header=false))
    winddf = DataFrame(CSV.File(winddatafile, skipto=2, header=false))

    println("compare")
    println(sum.(eachcol(baseffdf .- basesowfadf)))
    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    nturbines = length(basesowfadf.Column1)
    println("nturbs = ", nturbines)

    ndirections = length(winddf.d)

    # compute directional data 
    if normalize
        normalization_factor_sowfa = nturbines*maximum(maximum(eachrow(basesowfadf)))
        normalization_factor_bp = nturbines*maximum(maximum(eachrow(baseffdf)))
    else
        normalization_factor_sowfa = 1E6
        normalization_factor_bp = 1E6
    end

    println(normalization_factor_sowfa)
    println(normalization_factor_bp)

    basedirpowerff = sum.(eachcol(baseffdf))./normalization_factor_bp
    optdirpowerff = sum.(eachcol(optffdf))./normalization_factor_bp
    basedirpowersowfa = sum.(eachcol(basesowfadf))./normalization_factor_sowfa
    optdirpowersowfa = sum.(eachcol(optsowfadf))./normalization_factor_sowfa

    # create directional power bar charts
    fig, ax = plt.subplots(1, figsize=[6,4], subplot_kw=Dict("projection" => "polar"))

    if plottype == "power"
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , optdirpowersowfa, color=colors[3], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("SOWFA-optimized", (160, 65), color=colors[3])
        case=="low-ti" && ax.annotate("SOWFA-optimized", (160, 62), color=colors[3])

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , basedirpowersowfa, color=colors[2], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("SOWFA-base", (160, 55.5), color=colors[2])
        case=="low-ti" && ax.annotate("SOWFA-base", (160, 54.5), color=colors[2])
        
        # ax.bar(winddf.d .+ 5, optdirpowerff.*1E-6, label="BP-Opt", width=10, color=colors[4])
        ax.plot(winddf.d , optdirpowerff, color=colors[4], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("BP-optimized", (160, 60), color=colors[4])
        case=="low-ti" && ax.annotate("BP-optimized", (160, 56.5), color=colors[4])

        # ax.bar(winddf.d .+ 5, basedirpowerff.*1E-6, label="BP-Base", width=10, color=colors[1])
        ax.plot(winddf.d , basedirpowerff, color=colors[1], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("BP-base", (160, 53), color=colors[1])
        case=="low-ti" && ax.annotate("BP-base", (160, 50), color=colors[1])

        # format the figure
        ax.set(xticks=winddf.d, ylim=[40, 70], xlabel="Direction (degrees)", ylabel="Directional power (MW)")
        ax.legend(frameon=false,ncol=2)
    elseif plottype == "improvement"
        
        improvementsowfa = (optdirpowersowfa./basedirpowersowfa .- 1.0)*100
        improvementff = (optdirpowerff./basedirpowerff .- 1.0)*100
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , improvementsowfa, color=colors[3], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("SOWFA", (100, 9), color=colors[3])
        case=="low-ti" && ax.annotate("SOWFA", (100, 10), color=colors[3])

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , improvementff, color=colors[2], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("BP", (160, 6), color=colors[2])
        case=="low-ti" && ax.annotate("BP", (160, 9), color=colors[2])

        # format the figure
        ax.set(xticks=winddf.d, ylim=[0, 15], xlabel="Direction (degrees)", ylabel="Directional Improvment (%)")
        ax.legend(frameon=false,ncol=2)

    elseif plottype == "error"
        errorbase = (basedirpowerff./basedirpowersowfa .- 1.0)*100
        erroropt = (optdirpowerff./optdirpowersowfa .- 1.0)*100
        # ax.bar(winddf.d .- 5, optdirpowersowfa.*1E-6, label="SOWFA-Opt", width=10, color=colors[3])
        ax.plot(winddf.d , errorbase, color=colors[3], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("Base", (130, -4.5), color=colors[3])
        case=="low-ti" && ax.annotate("Base", (130, -0.5), color=colors[3])

        # ax.bar(winddf.d .- 5, basedirpowersowfa.*1E-6, label="SOWFA-Base", width=10, color=colors[2])
        ax.plot(winddf.d , erroropt, color=colors[2], marker="o", linestyle="--")
        case=="high-ti" && ax.annotate("Optimized", (130, -8), color=colors[2])
        case=="low-ti" && ax.annotate("Optimized", (130, -2.6), color=colors[2])

        # format the figure
        ax.set(xticks=winddf.d, ylim=[-10, 2], xlabel="Direction (degrees)", ylabel="Directional Error (%)")
        ax.legend(frameon=false,ncol=2)
    elseif plottype == "wakeloss" || plottype == "annualenergyloss"

        # find which turbines are not waked
        # baseffarray = convert(Array{Float64, (nturbines, ndirections)}, baseffdf)
        # optwcarray = convert(Array, optwcdf)
        ideal_power_base_ff = zeros(ndirections)
        ideal_power_opt_ff = zeros(ndirections)
        ideal_power_base_sowfa = zeros(ndirections)
        ideal_power_opt_sowfa = zeros(ndirections)

        # find ideal power
        for j = eachindex(winddf.d)
            nb = 0
            no = 0
            for i = 1:nturbines
                if basewcdf[i,j] == 0.0
                    ideal_power_base_ff[j] += baseffdf[i,j]./normalization_factor_bp 
                    ideal_power_base_sowfa[j] += basesowfadf[i,j]./normalization_factor_sowfa
                    nb += 1
                end
                if optwcdf[i,j] == 0.0
                    ideal_power_opt_ff[j] += optffdf[i,j]./normalization_factor_bp
                    ideal_power_opt_sowfa[j] += optsowfadf[i,j]./normalization_factor_sowfa
                    no += 1
                end
            end
            ideal_power_base_ff[j] /= nb
            ideal_power_base_sowfa[j] /= nb
            ideal_power_opt_ff[j] /= no
            ideal_power_opt_sowfa[j] /= no
        end
        
        ideal_power_base_ff .*= nturbines
        ideal_power_base_sowfa .*= nturbines
        ideal_power_opt_ff .*= nturbines
        ideal_power_opt_sowfa .*= nturbines

        # get wake loss 
        basewakelossff = 100.0.*(1.0 .- basedirpowerff./ideal_power_base_ff)
        basewakelosssowfa = 100.0.*(1.0 .- basedirpowersowfa./ideal_power_base_sowfa)
        optwakelossff = 100.0.*(1.0 .- optdirpowerff./ideal_power_opt_ff)
        optwakelosssowfa = 100.0.*(1.0 .- optdirpowersowfa./ideal_power_opt_sowfa)

        if plottype == "annualenergyloss"
            ts = winddf.p*365*24.0.*1E-3.*(1/100)
            basewakelossff .*= ts.*ideal_power_base_ff
            basewakelosssowfa .*= ts.*ideal_power_base_sowfa
            optwakelossff .*= ts.*ideal_power_opt_ff
            optwakelosssowfa .*= ts.*ideal_power_opt_sowfa
            ts = ts[6]
        end

        # make data a complete loop
        basewakelossff = [basewakelossff; basewakelossff[1]]
        basewakelosssowfa = [basewakelosssowfa; basewakelosssowfa[1]]
        optwakelossff = [optwakelossff; optwakelossff[1]]
        optwakelosssowfa = [optwakelosssowfa; optwakelosssowfa[1]]
        push!(winddf, first(winddf))

        if plottype == "annualenergyloss"
            fticks = 0:5:20
            units = "GWh"
        elseif plottype == "wakeloss"
            fticks = 0:5:25
            units = "%"
        end
        # plot
        d2r = pi/180.0
        ff.plotwindrose!(ax, winddf.d*d2r, optwakelosssowfa, plotcommand="plot", kwargs=(color=colors[3], marker="o", linestyle="--", label="SOWFA-optimized"))
        
        ff.plotwindrose!(ax, winddf.d*d2r, basewakelosssowfa, plotcommand="plot", kwargs=(color=colors[2], marker="o", linestyle="--", label="SOWFA-base"))
        
        ff.plotwindrose!(ax, winddf.d*d2r, optwakelossff, plotcommand="plot", kwargs=(color=colors[4], marker="o", linestyle="--", label="BP-optimized"))
        
        ff.plotwindrose!(ax, winddf.d*d2r, basewakelossff, plotcommand="plot", fticks=fticks, units=units, kwargs=(color=colors[1], marker="o", linestyle="--", label="BP-base"))
        
        ax.legend(frameon=false,ncol=1,loc=0,bbox_to_anchor=(1, 1))
        ax.set_axisbelow(false)
    end

    # remove upper and right bounding box
    # ax.spines["right"].set_visible(false)
    # ax.spines["top"].set_visible(false)

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

# create compound figure with all directional results 
function directional_comparison_figure_compound(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="directional-comparison")

    fig, ax = plt.subplots(3, 2, figsize=[8,9], sharex=true, sharey="row")

    # power
    directional_comparison_figure(colors, fontsize, ax=ax[1,1], case="high-ti", tuning="sowfa-nrel", plottype="power")
    ax[1,1].annotate("SOWFA-optimized", (160, 62), color=colors[3])
    ax[1,1].annotate("SOWFA-base", (150, 55.85), color=colors[2])
    ax[1,1].annotate("BP-optimized", (160, 60), color=colors[4])
    ax[1,1].annotate("BP-base", (160, 52.75), color=colors[1])
    directional_comparison_figure(colors, fontsize, ax=ax[1,2], case="low-ti", tuning="sowfa-nrel", plottype="power")
    ax[1,2].annotate("SOWFA-optimized", (160, 61), color=colors[3])
    ax[1,2].annotate("SOWFA-base", (150, 55), color=colors[2])
    ax[1,2].annotate("BP-optimized", (160, 57.5), color=colors[4])
    ax[1,2].annotate("BP-base", (80, 53), color=colors[1])
        
    # Error
    directional_comparison_figure(colors, fontsize, ax=ax[2,1], case="high-ti", tuning="sowfa-nrel", plottype="error")
    ax[2,1].annotate("Base", (130, -4.5), color=colors[3])
    ax[2,1].annotate("Optimized", (130, -6.5), color=colors[2])
    directional_comparison_figure(colors, fontsize, ax=ax[2,2], case="low-ti", tuning="sowfa-nrel", plottype="error")
    ax[2,2].annotate("Base", (130, -0.4), color=colors[3])
    ax[2,2].annotate("Optimized", (130, -2.7), color=colors[2])

    # Improvement
    directional_comparison_figure(colors, fontsize, ax=ax[3,1], case="high-ti", tuning="sowfa-nrel", plottype="improvement")
    ax[3,1].annotate("SOWFA", (91, 9), color=colors[3])
    ax[3,1].annotate("BP", (160, 6.8), color=colors[2])
    directional_comparison_figure(colors, fontsize, ax=ax[3,2], case="low-ti", tuning="sowfa-nrel", plottype="improvement")
    ax[3,2].annotate("SOWFA", (83, 10), color=colors[3])
    ax[3,2].annotate("BP", (160, 8.8), color=colors[2])

    # add x label 
    ax[3,1].set(xlabel="Direction (degrees)", xticks=10:60:360)
    ax[3,2].set(xlabel="Direction (degrees)", xticks=10:60:360)

    # add y label 
    ax[1,1].set(ylabel="Power (MW)", xticks=10:60:360, yticks=50:2:66, ylim=[50,66])
    ax[1,2].set(xticks=10:60:360, yticks=50:2:66)
    ax[2,1].set(ylabel="Power Error (%)", xticks=10:60:360, yticks=-10:2:2, ylim=[-10,2])
    ax[2,2].set(xticks=10:60:360, yticks=-10:2:2)
    ax[3,1].set(ylabel="Power Improvement (%)", xticks=10:60:360, yticks=0:2:16, ylim=[0,16])
    ax[3,1].set(xticks=10:60:360, yticks=0:2:16)

    # add column titles
    ax[1,1].set(title="High-TI")
    ax[1,2].set(title="Low-TI")

    # add subfigure labels
    i=97
    for axi in [ax[1,:]; ax[2,:]; ax[3,:]]
        axi.annotate("($(Char(i)))", xy=(0.05, 0.05), xycoords="axes fraction")
        i += 1
    end
    plt.tight_layout()

    # save figure
    if showfigs
        plt.show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*"power"*".pdf", transparent=true)
    end

    fig, ax = plt.subplots(3, 2, figsize=[8,9], sharex="col", sharey="row")
    # Wake Loss
    directional_comparison_figure(colors, fontsize, ax=ax[1,1], case="high-ti", tuning="sowfa-nrel", plottype="wakeloss")
    ax[1,1].annotate("SOWFA-optimized", (21, 5.5), color=colors[3])
    ax[1,1].annotate("SOWFA-base", (155, 13.5), color=colors[2])
    ax[1,1].annotate("BP-optimized", (49.6, 9.9), color=colors[4])
    ax[1,1].annotate("BP-base", (97, 18), color=colors[1])
    directional_comparison_figure(colors, fontsize, ax=ax[1,2], case="low-ti", tuning="sowfa-nrel", plottype="wakeloss")
    ax[1,2].annotate("SOWFA-optimized", (160, 9.2), color=colors[3])
    ax[1,2].annotate("SOWFA-base", (170, 18), color=colors[2])
    ax[1,2].annotate("BP-optimized", (113, 16), color=colors[4])
    ax[1,2].annotate("BP-base", (95, 22), color=colors[1])

    # Annual Directional Energy Loss
    directional_comparison_figure(colors, fontsize, ax=ax[2,1], case="high-ti", tuning="sowfa-nrel", plottype="annualenergyloss")
    directional_comparison_figure(colors, fontsize, ax=ax[2,2], case="low-ti", tuning="sowfa-nrel",  plottype="annualenergyloss")

    # Change in Directional Energy Loss
    directional_comparison_figure(colors, fontsize, ax=ax[3,1], case="high-ti", tuning="sowfa-nrel", plottype="annualenergylossincrease")
    ax[3,1].annotate("SOWFA", (130, 5), color=colors[3])
    ax[3,1].annotate("BP", (212, 5), color=colors[2])
    directional_comparison_figure(colors, fontsize, ax=ax[3,2], case="low-ti", tuning="sowfa-nrel",  plottype="annualenergylossincrease")
    ax[3,2].annotate("SOWFA", (122, 5), color=colors[3])
    ax[3,2].annotate("BP", (199, 5), color=colors[2])

    # add x labels and ticks
    ax[1,1].set(xticks=10:60:360)
    ax[1,2].set(xticks=10:60:360)
    ax[2,1].set(xticks=10:60:360)
    ax[2,2].set(xticks=10:60:360)
    ax[3,1].set(xlabel="Direction (degrees)", xticks=10:60:360)
    ax[3,2].set(xlabel="Direction (degrees)", xticks=10:60:360)

    # add y labels and ticks
    ax[1,1].set(ylabel="Wake Loss (%)", yticks=0:5:30, ylim=[0, 30])
    ax[1,2].set(yticks=0:5:30, ylim=[0, 30])
    ax[2,1].set(ylabel="Annual Energy Loss (GW h)", yticks=0:5:20, ylim=[0, 20])
    ax[2,2].set(yticks=0:5:20, ylim=[0, 20])
    ax[3,1].set(ylabel="Annual Energy Improvement (GW h)", yticks=0:2:12, ylim=[0, 12])
    ax[3,2].set(yticks=0:2:12, ylim=[0, 12])

    # add legend
    ax[2,1].legend(frameon=false)
    ax[2,2].legend(frameon=false)

    # add column titles
    ax[1,1].set(title="High-TI")
    ax[1,2].set(title="Low-TI")

    # add subfigure labels
    i=97
    for axi in [ax[1,:]; ax[2,:]; ax[3,:]]
        axi.annotate("($(Char(i)))", xy=(0.05, 0.05), xycoords="axes fraction")
        i += 1
    end

    plt.tight_layout()

    # save figure
    if showfigs
        plt.show()
    end
    if savefigs
        plt.savefig(image_directory*image_name*"wakeloss"*".pdf", transparent=true)
    end

end

function turbine_comparison_figures(colors, fontsize; ax1=nothing, ax2=nothing, showfigs=false, savefigs=false, image_directory="images/", image_name="turbine-comparison", case="low-ti", tuning="alldirections")
    
    function plot_turbine_heatmap(data, winddirections, vmin, vmax; ax=nothing, wake_count=nothing, drawcbar=true)

        # number of turbines in the farm
        nturbines = 38

        # number of wind directions 
        ndirections = length(winddirections)

        # intialize figure
        if ax === nothing
            fig, ax = plt.subplots(figsize=(8,4))
        end

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
                cbarlabel="Turbine power error as percent of maximum SOWFA turbine power", cbar_kw=d,
                use_cbar=drawcbar)

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

    if ax1 !== nothing
        drawcbar2 = false
    else
        drawcbar2 = true
    end

    if case == "high-ti"
        layoutid = 83
        n = 4
    elseif case == "low-ti"
        layoutid = 252
        n = 2
    end

    # load wake data 
    wakecountfile = "image-data/wakes/turbine-wakes-opt$n-ff-$case-$tuning-layout$layoutid.txt"
    
    # flowfarm base data 
    basepowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout1-base.txt"

    # sowfa base data
    basepowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les.txt"

    # flowfarm opt data 
    optpowerfileff = "image-data/power/FLOWFarm/turbine-power-ff-100pts-$case-$tuning-layout$layoutid-opt$n.txt"
    
    # sowfa opt data 
    optpowerfilesowfa = "image-data/power/SOWFA/turbine-power-$case-les-opt$n.txt"

    # windrose data 
    winddatafile = "../src/inputfiles/wind/windrose_nantucket_12dir.txt"

    # read files to dataframes
    wake_count = transpose(readdlm(wakecountfile, ',', skipstart=1, header=false))
    baseff = transpose(readdlm(basepowerfileff, ',', skipstart=1, header=false))
        # DataFrame(CSV.File(basepowerfileff, skipto=2, header=false))
    optff = transpose(readdlm(optpowerfileff, ',', skipstart=1, header=false))
    basesowfa = transpose(readdlm(basepowerfilesowfa, skipstart=1, header=false))
    # DataFrame(CSV.File(basepowerfilesowfa, skipto=2, header=false))
    optsowfa = transpose(readdlm(optpowerfilesowfa, skipstart=1, header=false))
    winddf = DataFrame(CSV.File(winddatafile, skipto=2, header=false))

    # name wind data columns 
    rename!(winddf,:Column1 => :d,:Column2 => :s,:Column3 => :p)

    # set vmin and vmax 
    vmax = 28
    vmin = -vmax

    # calculate turbine errors for base case
    turberror = errors(basesowfa, baseff, method="normbyfirst")
    ers = (sum(basesowfa[i,:] for i=1:12 ) .- sum(baseff[j,:] for j=1:12))
    # ers = (basesowfa .- baseff)
    println("base turbine error sum $case $tuning: $(sum(turberror))")
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
    plot_turbine_heatmap(data, winddf.d, vmin, vmax; ax=ax1)

    if savefigs
        if !drawcbar2
            plt.savefig(image_directory*image_name*"-"*case*"-"*tuning*".pdf", transparent=true)
        end
    end

    # calculate turbine errors for opt case
    turberror = errors(optsowfa, optff, method="normbyfirst")
    println("opt turbine error sum $case $tuning: $(sum(turberror))")
    data = convert.(Int64, round.(turberror.*100, digits=0))

    # plot error on heatmap
    plot_turbine_heatmap(data, winddf.d, vmin, vmax, ax=ax2, drawcbar=drawcbar2)

    if savefigs
        if !drawcbar2
            plt.savefig(image_directory*image_name*"-"*case*"-"*tuning*"-opt.pdf", transparent=true)
        end
    end

    # save figure
    if showfigs
        plt.show()
    end

end

function turbine_comparison_figures_compound(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="turbine-comparison", case="low-ti", tuning="alldirections")
    
    # generate figure
    fig, ax = plt.subplots(2,1,figsize=(8,8))
    
    # populate figure
    turbine_comparison_figures(colors, fontsize, ax1=ax[1], ax2=ax[2], savefigs=savefigs, showfigs=showfigs, case=case, tuning="sowfa-nrel")
    
    # format figure
    ax[1].set(ylabel="Wind Direction (degrees)", xlabel="Turbine Index")
    ax[1].xaxis.set_label_position("top") 
    if case == "high-ti"
        ax[1].set_title("(a) High-TI Base")
    elseif case == "low-ti"
        ax[1].set_title("(a) Low-TI Base")
    end
    
    ax[2].set(ylabel="Wind Direction (degrees)", xlabel="Turbine Index")
    ax[2].xaxis.set_label_position("top") 
    if case == "high-ti"
        ax[2].set_title("(b) High-TI Optimized")
    elseif case == "low-ti"
        ax[2].set_title("(b) Low-TI Optimized")
    end

    plt.tight_layout()

    if savefigs
        plt.savefig(image_directory*image_name*"-"*case*"-"*tuning*"-compound.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function horns_rev_rows_verification_figure(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="horns-rev-rows")


    # this function is called to populate each subfigure
    function each_subfig!(ax, nsamplepoints, colors, fontsize)
        # load FLOWFarm and Niayifar data
        datafile = "image-data/verification/horns-rev-rows-$nsamplepoints-sample-points.txt"
        data = readdlm(datafile, skipstart=1)
        rows = data[:, 1]
        normalized_power_les_niayifar = data[:,2]
        normalized_power_model_niayifar = data[:, 3]
        normalized_power_averaged_ff_no_ti = data[:,4]
        normalized_power_averaged_ff_ti = data[:, 5]

        markersize = 15
        ax.scatter(rows, normalized_power_les_niayifar, c=colors[1], label="Niayifar 2016 LES", marker="o", s=markersize)
        ax.scatter(rows, normalized_power_model_niayifar, c=colors[4], label="Niayifar 2016 model", marker="*", s=markersize)
        ax.scatter(rows, normalized_power_averaged_ff_no_ti, edgecolors=colors[2], label="BP w/o local TI", marker="^", facecolors="none", s=markersize)
        ax.scatter(rows, normalized_power_averaged_ff_ti, edgecolors=colors[3], label="BP w/local TI", marker="v", facecolors="none", s=markersize)
        
        # format the figure
        ax.set(xlabel="Row", ylabel="Normalized power", ylim=[0,1.1], xticks=Int.(rows))
        ax.legend(frameon=false,ncol=1)

        # remove upper and right bounding box
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        # calculate the appropriate x location for the given figure 
        xtxtloc = (maximum(rows) + minimum(rows))/2
        println(xtxtloc)
        return xtxtloc
    end

    fig, ax = plt.subplots(1,2,figsize=(10,3))

    xtxtloc100 = each_subfig!(ax[1], 100, colors, fontsize)
    xtxtloc1 = each_subfig!(ax[2], 1, colors, fontsize)

    ax[1].text(xtxtloc100,-0.3,"(a)",horizontalalignment="center",fontsize=fontsize*1.25)
    ax[2].text(xtxtloc1,-0.3,"(b)",horizontalalignment="center",fontsize=fontsize*1.25)

    plt.subplots_adjust(top=0.99,bottom=0.22,right=0.92,left=0.08,wspace=5.4)

    plt.tight_layout()

    if savefigs
        plt.savefig(image_directory*image_name*".pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function horns_rev_direction_verification_figure(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="horns-rev-directions")
    
    # this function is called to populate each subfigure
    function each_subfig!(ax, nsamplepoints, colors)
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

        ax.plot(lesdirections, normalized_power_les_niayifar, c=colors[1], label="Niayifar 2016 LES", marker="o", linestyle="none", markersize=3)
        ax.plot(modeldirections, normalized_power_model_niayifar, "-", c=colors[4], label="Niayifar 2016 model")
        ax.plot(modeldirections, normalized_power_averaged_ff_no_ti, "--", c=colors[2], label="BP w/o local TI")
        ax.plot(modeldirections, normalized_power_averaged_ff_ti, ":", c=colors[3], label="BP w/local TI")
        
        # format the figure
        ax.set(xlabel="Direction (degrees)", ylabel="Normalized power", ylim=[0,1])
        ax.legend(frameon=false, ncol=2)

        # remove upper and right bounding box
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        # calculate the appropriate x location for the given figure 
        xtxtloc = (maximum(modeldirections) + minimum(modeldirections))/2
        println(xtxtloc)
        return xtxtloc
    end

    fig, ax = plt.subplots(1,2, figsize=(10,3))
    
    xtxtloc100 = each_subfig!(ax[1], 100, colors)
    xtxtloc1 = each_subfig!(ax[2], 1, colors)

    ax[1].text(xtxtloc100,-0.3,"(a)",horizontalalignment="center", fontsize=fontsize*1.25)
    ax[2].text(xtxtloc1,-0.3,"(b)",horizontalalignment="center", fontsize=fontsize*1.25)

    plt.subplots_adjust(top=0.99,bottom=0.22,right=0.92,left=0.08,wspace=0.8)

    plt.tight_layout()

    if savefigs
        plt.savefig(image_directory*image_name*".pdf", transparent=true)
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


    ax1.plot(x1,y1,"o",color=color,markersize=4,label="Sampling points")
    ax2.plot(x2,y2,"o",color=color,markersize=4)
    ax1.legend(fontsize=fontsize)

    # ax1.set_xlabel(r"Horizontal Distance From Hub, $\Delta y/d$",fontsize=fontsize)
    ax1.set_xlabel("Horizontal distance from hub, "*L"\Delta y/d",fontsize=fontsize)
    ax1.set_ylabel("Vertical distance from hub, "*L"\Delta z/d",fontsize=fontsize)
    ax2.set_xlabel("Horizontal distance from hub, "*L"\Delta y/d",fontsize=fontsize)
    ax2.set_ylabel("Vertical distance from hub, "*L"\Delta z/d",fontsize=fontsize)

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

function turbine_layouts(colors ;ax=nothing, rotor_diameter=126.4,les_side=5000,
                                    c1="C0",c2="C1",color="C0",fontsize=10,
                                    case="low-ti", tuning="sowfa-nrel", showfigs=false, savefigs=false,
                                    iter="base", layoutid=1, n=1, gen="angle-each-circle",
                                    compound=true, annotate=false)

    les_radius = les_side/2.0
    boundary_radius = les_radius - 500.0 - rotor_diameter/2.0

    # load data
    if iter == "base"
        datafile = "../src/inputfiles/farms/layout_38turb_round.txt"
        data = readdlm(datafile, skipstart=1).*rotor_diameter .+ (les_radius - boundary_radius)
    elseif iter == "start"
        # datafile = "image-data/layouts/opt/optresultsmilestone.csv"
        datafile = "../src/inputfiles/farms/startinglayouts/$gen/nTurbs38_spacing5.0_layout_$layoutid.txt"
        data = readdlm(datafile, ',', skipstart=1).*rotor_diameter .+ (les_radius)
    elseif iter == "opt"
        # datafile = "image-data/layouts/opt/optresultsmilestone.csv"
        datafile = "image-data/layouts/opt/$case-$tuning-opt$n-layout$layoutid-aec-wec.csv"
        data = readdlm(datafile, ',', skipstart=1) .+ les_radius
    else
        ErrorException("Case not available")
    end

    xlocs = data[:,1]
    ylocs = data[:,2]

    if ax === nothing
        plt.figure(figsize=(5,3))
        ax = plt.gca()
    end
    
    if !compound
        ax.axes.xaxis.set_visible(false)
        ax.axes.yaxis.set_visible(false)

        ax.spines["top"].set_visible(false)
        ax.spines["bottom"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["right"].set_visible(false)
        plt.subplots_adjust(left=0.15,bottom=0.12,top=0.99,right=0.6)
    end
    if annotate == true
        ax.annotate("Farm boundary", (les_radius/2, les_radius-boundary_radius-rotor_diameter*3.0), color=colors[2])
        ax.annotate("LES domain", (0, les_side+rotor_diameter/2.0), color=colors[4])

    end

    # xaxis.set_visible(False)

    R = rotor_diameter/2
    cx = 2500
    cy = 2500
    lw = 0.75
    plot_circle(cx,cy,boundary_radius,colors[2],ax,linestyle="--",linewidth=lw,label="Farm boundary")
    
    les_x = [0,les_side,les_side,0,0]
    les_y = [0,0,les_side,les_side,0]
    ax.plot(les_x,les_y,"-",color=colors[4],linewidth=lw,label="LES domain")
    
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
    
   

    plt.tight_layout()

    if savefigs
        plt.savefig("images/"*case*"-layout.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function turbine_layouts_compound(colors ;rotor_diameter=126.4,les_side=5000,
    c1="C0",c2="C1",color="C0",fontsize=10,
    case="low-ti-opt", showfigs=false, savefigs=false)

    fig, ax = plt.subplots(2,3, figsize=(8,6), gridspec_kw=Dict("wspace" => 0.05, "hspace"=>0.0))

    # high ti
    turbine_layouts(colors ;ax=ax[1,1], fontsize=fontsize, iter="base", layoutid=1, case="high-ti", annotate=true)
    turbine_layouts(colors ;ax=ax[1,2], fontsize=fontsize, iter="start", layoutid=83, n=4, case="high-ti")
    turbine_layouts(colors ;ax=ax[1,3], fontsize=fontsize, iter="opt", layoutid=83, n=4, case="high-ti")
    ax[1,1].set(ylabel="High-TI")
    
    # low ti
    turbine_layouts(colors ;ax=ax[2,1], fontsize=fontsize, iter="base", layoutid=1, case="low-ti")
    turbine_layouts(colors ;ax=ax[2,2], fontsize=fontsize, iter="start", layoutid=252, n=2, case="low-ti")
    turbine_layouts(colors ;ax=ax[2,3], fontsize=fontsize, iter="opt", layoutid=252, n=2, case="low-ti")
    ax[2,1].set(ylabel="Low-TI")

    # add overhead labels 
    ax[1,1].set(title="Base")
    ax[1,2].set(title="Start")
    ax[1,3].set(title="Optimized")

    # remove spines
    for axi in [ax[1,:]; ax[2,:]]
        axi.axes.xaxis.set_visible(false)
        # axi.axes.yaxis.set_visible(false)
        axi.tick_params(left=false, labelleft=false)
        #remove background patch (only needed for non-white background)
        axi.patch.set_visible(false)

        axi.spines["top"].set_visible(false)
        axi.spines["bottom"].set_visible(false)
        axi.spines["left"].set_visible(false)
        axi.spines["right"].set_visible(false)
    end

    # add subfigure labels
    i=97
    for axi in [ax[1,:]; ax[2,:]]
        axi.annotate("($(Char(i)))", xy=(0.1, 0.1), xycoords="axes fraction")
        i += 1
    end

    plt.tight_layout()

    if savefigs
        plt.savefig("images/layouts-compound.pdf", transparent=true)
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

function vertical_slice(colors; savefigs=false, showfigs=false, fontsize=10)

    # for the high-ti case 

    # set input values 
    diam = 126.4
    hub_height = 90.0

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
    levels = 0:0.5:10

    # get colormap 
    cmap = custum_color_map(idx=[3,2])
    
    # generate contour plot
    fig, ax = plt.subplots(2, figsize=(6,4))
    
    ffvelocities[ffvelocities .< 0] .= 0
    cs = ax[1].contourf(xg./diam, zg./diam, ffvelocities, levels, cmap=cmap)
    ax[2].contourf(xg./diam, zg./diam, ffvelocities, levels, cmap=cmap)
    
    # cover up interpolated region for non-interpolated figure
    # r = plt.matplotlib.patches.Rectangle((0,0),424.39110240363897/diam,2, color="w")
    r = plt.matplotlib.patches.Rectangle((0,0),305.94115534342194/diam,2, color="w")
    ax[1].add_patch(r)

    position=fig.add_axes([0.8,0.2,0.03,0.725])
    cbar = fig.colorbar(cs, cax=position, ax=[ax[1],ax[2]], label=L"Wind speed (m s$^{-1}$)", use_gridspec=true)
    
    # add turbine 
    radius = 0.5
    thub = hub_height/diam
    hd = 0.1
    chord = 0.2
    nacellewidth = 3*chord
    nacelleheigt = hd
    towerdtop = 3.87/diam 
    towerdbot = 6.0/diam
    overhang = 5.0/diam

    ff.add_turbine!(ax[1], view="side", hubdiameter=0.1, hubheight=hub_height/diam, radius=0.5, chord=0.2, nacellewidth=0.6, nacelleheight=0.1, towerbottomdiam=6.0/diam, towertopdiam=3.87/diam, overhang=5.0/diam, s=5)

    ff.add_turbine!(ax[2], view="side", hubdiameter=0.1, hubheight=hub_height/diam, radius=0.5, chord=0.2, nacellewidth=0.6, nacelleheight=0.1, towerbottomdiam=6.0/diam, towertopdiam=3.87/diam, overhang=5.0/diam, s=5)

    # add labels
    # ax[1].set_title("(a)", fontsize=fontsize, y=0)#text(0.5,-0, "(a)", horizontalalignment="center", fontsize=fontsize, transform=ax[1].transAxes)
    # plt.annotate( "(b)", (0.5,-1), xycoords=fig.transFigure, horizontalalignment="center", fontsize=fontsize)
    figt = plt.gcf()
    xtxtloc = (xmin + (xmax-xmin)/2.0)/diam
    
    ax[1].text(xtxtloc,-1.,"(a)",horizontalalignment="center")
    ax[2].text(xtxtloc,-1.,"(b)",horizontalalignment="center")

    plt.subplots_adjust(top=1,bottom=0.1,right=0.92,left=0.12,hspace=0.5)

    # format figure 
    ax[1].set(xticks=0:4:20, yticks=0:2)
    ax[2].set(xticks=0:4:20, yticks=0:2)
    ax[1].set(xlabel=L"$x/d_0$", ylabel=L"$z/d_0$")
    ax[2].set(xlabel=L"$x/d_0$", ylabel=L"$z/d_0$")

    plt.subplots_adjust(hspace=0.8)

    plt.tight_layout(pad=0.5, rect=(0, 0, 0.75, 1))

    if savefigs
        plt.savefig("images/vertical-slice.pdf", transparent=true)
    end

    # save figure
    if showfigs
        plt.show()
    end
end

function plot_results_distribution(colors; savefigs=false, showfigs=false, fontsize=10)

    # load data 
    high_ti_data = DataFrame(CSV.File("./image-data/power/multi-start-distributions/results-distribution-high-ti-4.csv", skipto=2))
    low_ti_data = DataFrame(CSV.File("./image-data/power/multi-start-distributions/results-distribution-low-ti-1.csv", skipto=2))
    
    # generate figure
    fig, ax = plt.subplots(1,2, sharey=true, sharex=true, figsize=(8,3))

    # plot histogram 
    bins = 400:1:500
    (n, binedges, patches) = ax[1].hist(high_ti_data[!, :aepib], label="Start AEP", color=colors[1],bins=bins)
    bx = high_ti_data[1, :aepib]
    by = Int(n[findall(<=(high_ti_data[1, :aepib]), binedges)[end]])
    ax[1].hist(ones(by)*bx, label="Base Start Bin", color=colors[1], edgecolor=colors[2], bins=bins)

    ax[1].hist(high_ti_data[!, :aepfb], label="Optimized AEP", color=colors[3],bins=bins)
    
    (n, binedges, patches) = ax[2].hist(low_ti_data[!, :aepib]*1E-1, label="Start AEP", color=colors[1],bins=bins)
    bx = low_ti_data[1, :aepib]*1E-1
    by = Int(n[findall(<=(low_ti_data[1, :aepib]*1E-1), binedges)[end]])
    ax[2].hist(ones(by)*bx, label="Base Start Bin", color=colors[1], bins=bins, edgecolor=colors[2])

    ax[2].hist(low_ti_data[!, :aepfb]*1E-1, label="Optimized AEP", color=colors[3],bins=bins)
    
    # format histogram
    ax[1].set(xlim=[400,500], ylim=[0,80], xlabel="AEP (GW h)", ylabel="Count", title="(a) High-TI")
    ax[2].set(xlim=[400,500], ylim=[0,80], xlabel="AEP (GW h)", title="(b) Low-TI")
    ax[1].legend(frameon=false)
    ax[2].legend(frameon=false, loc="upper left")

    # removes spines
    ax[1].spines["right"].set_visible(false)
    ax[1].spines["top"].set_visible(false)
    ax[1].spines["bottom"].set_visible(true)
    ax[1].spines["left"].set_visible(true)

    ax[2].spines["right"].set_visible(false)
    ax[2].spines["top"].set_visible(false)
    ax[2].spines["bottom"].set_visible(true)
    ax[2].spines["left"].set_visible(false)
    ax[2].tick_params(
    axis="y",          # changes apply to the x-axis
    which="both",      # both major and minor ticks are affected
    left=false,      # ticks along the bottom edge are off
    right=false,         # ticks along the top edge are off
    labelleft=false) # labels along the bottom edge are off

    plt.tight_layout()

    # save figure
    if savefigs
        plt.savefig("images/aep-distributions-compound.pdf", transparent=true)
    end

    # show figure
    if showfigs
        plt.show()
    end
end

function make_images()

    fontsize = 8
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]

    rcParams = PyPlot.matplotlib.rcParams
    rcParams["font.size"] = fontsize
    rcParams["lines.markersize"] = 1
    rcParams["axes.prop_cycle"] = colors

    savefigs = true 
    showfigs = true

    # wind_shear_tuning_compound(colors, fontsize, showfigs=showfigs, savefigs=savefigs, image_directory="images/", image_name="windshear-compound")

    # layout(colors, fontsize)

    # opt_comparison_table(case="high-ti", tuning="sowfa-nrel")
    # opt_comparison_table(case="low-ti", tuning="sowfa-nrel")
    
    # directional_comparison_figure_compound(colors, fontsize, showfigs=showfigs, savefigs=savefigs)

    # turbine_comparison_figures_compound(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel")
    # turbine_comparison_figures_compound(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel")
    
    plot_results_distribution(colors; savefigs=savefigs, showfigs=showfigs, fontsize=fontsize)

    # horns_rev_rows_verification_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs)
    # horns_rev_direction_verification_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs)

    # turbine_layouts_compound(colors, case="low-ti", showfigs=showfigs, savefigs=savefigs)

    # vertical_slice(colors, savefigs=savefigs, showfigs=showfigs)

    # windrose(d1,f1,d2,f2;color="C0",alpha=0.5,fontsize=8,filename="nosave")
    
end


    ############# obsolete/not used #############
    # wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti")
    # wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti")
    
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="power")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel", normalize=false, plottype="power")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="improvement")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel", normalize=false, plottype="improvement")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="error")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel", normalize=false, plottype="error")
    
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="wakeloss")
    # directional_comparison_figure(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel", normalize=false, plottype="wakeloss")
    
    # directional_comparison_figure_polar(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel", normalize=false, plottype="annualenergyloss")
    # directional_comparison_figure_polar(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel", normalize=false, plottype="annualenergyloss")
    

    # turbine_layouts(colors, case="low-ti-opt", showfigs=showfigs, savefigs=savefigs)
    # turbine_layouts(colors, case="low-ti", showfigs=showfigs, savefigs=savefigs)

    # turbine_layouts(colors, case="high-ti-opt", showfigs=showfigs, savefigs=savefigs)
    # turbine_layouts(colors, case="low-ti", showfigs=showfigs, savefigs=savefigs)

    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="low-ti", tuning="sowfa-nrel")
    # turbine_comparison_figures(colors, fontsize, savefigs=savefigs, showfigs=showfigs, case="high-ti", tuning="sowfa-nrel")
    