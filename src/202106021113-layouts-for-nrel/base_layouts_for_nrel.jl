using FLOWFarm; const ff=FLOWFarm
using DataFrames
using CSV
using DelimitedFiles
using PyPlot

function center_turbine_locations!(turbine_x, turbine_y, center)
    # turbine_x .-= turbine_x[1]
    # turbine_y .-= turbine_y[1]
    turbine_x .+= center[1]
    turbine_y .+= center[2]
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

function print_optimized_layouts(;df=nothing, outfile="optlayout", diam=126.4)
    if df === nothing
        include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs-low-ti.jl")
        df = DataFrame(CSV.File("../202107012007-optimization/optresultseasy.csv"))
        outfile = "round_38_turbines_opt.xlsx"
    end
    # println(df.x[end][1])
    turbine_x = df.x
    turbine_y = df.y
    println(typeof(df.x))
    rotor_diameter = zeros(length(turbine_x)) .+ diam
    center = [2500.0, 2500.0]
    center_turbine_locations!(turbine_x, turbine_y, center)
    ff.print_layout_in_cartesian_frame_excel(turbine_x, turbine_y, rotor_diameter, winddirections, outfile, center=center, round_locations=true, plot_layouts=true)
    
end