using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles

function generate_layouts(nlayouts; startingindex=1, method="individual")

    
    # specify directory 
    output_directory = "../inputfiles/farms/startinglayouts/$(method)2/"
    println(output_directory)
    show = true 
    save = false

    # specify parameters 
    rotor_diameter = 126.4
    farm_center = [0.0, 0.0]
    farm_diameter = 4000.0
    base_spacing = 5.0 
    min_spacing = 5.0

    # get layouts
    ff.generate_round_layouts(nlayouts, rotor_diameter, farm_center=farm_center, method=method,
                            farm_diameter=farm_diameter, base_spacing=base_spacing, min_spacing=min_spacing,
                            output_directory=output_directory, show=show, save_layouts=save, startingindex=startingindex)
end

function check_base_equality()

    orig = readdlm("../inputfiles/farms/layout_38turb_round.txt", skipstart=1)
    new = readdlm("../inputfiles/farms/startinglayouts/nTurbs38_spacing5.0_layout_1.txt", ',', skipstart=1)
    
    # println(new)
    orig[:,1] .-= orig[1,1]
    orig[:,2] .-= orig[1,2]

    # println(orig)
    println(maximum(new .- orig))

end