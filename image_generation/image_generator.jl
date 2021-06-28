import PyPlot; const plt = PyPlot
using DataFrames 
using CSV
using PyPlot, Color

function custum_color_map()
    colors = [colorant"#BDB8AD", colorant"#85C0F9", colorant"#0F2080", colorant"#F5793A", colorant"#A95AA1", colorant"#382119"]
    # @pyimport matplotlib.colors as matcolors
    # cmap = matcolors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],"A")

    return ColorMap("BlueGrayOrange", [colors[3],colors[1],colors[4]])
end

function heatmap(data, row_labels, col_labels; ax=nothing, cbar_kw=Dict(), cbarlabel="", use_cbar=true, labelpixels=true, fontsize=10, vcolor="w", edgecolor="w")
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

    # Plot the heatmap
    im = ax.pcolormesh(data, edgecolor=edgecolor, cmap=cbar_kw[:cmap], vmin=minimum(cbar_kw[:ticks]), vmax=maximum(cbar_kw[:ticks]))

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
                        ha="center",va="center",
                        size=fontsize,color=vcolor)
            end
        end
    end

    return im, cbar
end

function turbine_error(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="turbine-error-", case="lowti", layout="opt")
    image_name *= case*"-"*layout*".pdf"
    
end

function wind_shear_tuning(colors, fontsize; showfigs=false, savefigs=false, image_directory="images/", image_name="windshear.pdf")

    # load data 
    df_model = DataFrame(CSV.File("./image_data/wind_shear/wind_shear_tuned.txt", datarow=2, header=false))
    df_les = DataFrame(CSV.File("./image_data/wind_shear/wind_shear_les.txt", datarow=2, header=false))

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
    ax.annotate("LES", (7.75, 125), color=colors[2], alpha=1.0, size=fontsize)

    # plot model 
    ax.plot(df_model.s, df_model.h, color=colors[3])
    ax.annotate("Curve Fit", (8.15, 95), color=colors[3], alpha=1.0, size=fontsize)

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
        plt.savefig(image_directory*image_name, transparent=true)
    end
end
function generate_images_for_publication()
    fontsize = 18
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]
    savefigs = true 
    showfigs = true

    wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs)

end