import PyPlot; const plt = PyPlot
using DataFrames 
using CSV

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