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
    # uref=7.84 
    # shear=0.086 
    # ground=0.0


    # plot les data 
    p = plt.scatter(df_les.s, df_les.h, color=colors[2])

    # plot model 
    plt.plot(df_model.s, df_model.h, color=colors[3])

    # refine figure
    plt.xlabel("Wind Speed (m/s)")
    plt.ylabel("Height (m)")

    # save figure
    show(p)

end
function generate_images_for_publication()
    fontsize = 13
    colors = ["#BDB8AD", "#85C0F9", "#0F2080", "#F5793A", "#A95AA1", "#382119"]
    savefigs = true 
    showfigs = true

    wind_shear_tuning(colors, fontsize, savefigs=savefigs, showfigs=showfigs)

end