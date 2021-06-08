using FLOWFarm; const ff=FLOWFarm
using DataFrames
using Plots

function center_turbine_locations!(turbine_x, turbine_y, center)
    turbine_x .-= turbine_x[1]
    turbine_y .-= turbine_y[1]
end

function print_base_layouts()
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")
    center = [2500.0, 2500.0]
    # center_turbine_locations!(turbine_x, turbine_y, center)
    turbine_x .+= center[1]
    turbine_y .+= center[2]
    ff.print_layout_in_cartesian_frame_excel(turbine_x, turbine_y, winddirections, center=center, round_locations=true, plot_layouts=true)
    # p = plot()
    # ff.plotlayout!(p, turbine_x, turbine_y, rotor_diameter)
    # ff.plotlayout!(p, [turbine_x[1], turbine_x[20]], [turbine_y[1], turbine_y[20]], rotor_diameter, fillcolor=:red)
    # display(p)
end