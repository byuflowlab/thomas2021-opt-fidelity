using Base: windowserror
using CSV
include("../202105111413-SOWFA-comparison/202105111413-sowfa-round-farm-comparison.jl")

# compare turbine errors PlantEnergy and FLOWFarm 

function plantenergy_vs_flowfarm(;nturbines=38)
    # load plant energy data 
    # turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data(journal=false)
    turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019_new = get_data(journal=true)

    # load and run FLOWFarm setup 
    turbine_powers_by_direction_ff = run_flow_farm(use_local_ti=true, nsamplepoints=100, alpha=1.0, verbose=true, shearfirst=false)

    # calculate errors 
    difference_turbines_pe = turbine_powers_by_direction_sowfa - turbine_powers_by_direction_thomas2019_new

    # println(typeof(difference_turbines))

    # plant energy vs sowfa 
    df = DataFrame(SOWFA=turbine_powers_by_direction_sowfa[1,:], PlantEnergy=turbine_powers_by_direction_thomas2019_new[1,:], Diff=difference_turbines_pe[1,:])
    # df = DataFrame(PlantEnergy=turbine_powers_by_direction_thomas2019_new[1,:], FLOWFarm=turbine_powers_by_direction_thomas2019_new[1,:], Diff=difference_turbines_pe[1,:])
    println(describe(df))

    # flowfarm vs sowfa differences
    difference_turbines_ff = turbine_powers_by_direction_sowfa - turbine_powers_by_direction_ff
    dfff = DataFrame(SOWFA=turbine_powers_by_direction_sowfa[1,:], FLOWFarm=turbine_powers_by_direction_ff[1,:], Diff=difference_turbines_ff[1,:])
    println(describe(dfff))

    # flowfarm vs plant energy differences
    difference_turbines_ffvspe = turbine_powers_by_direction_thomas2019_new - turbine_powers_by_direction_ff
    dfffvspe = DataFrame(PlantEnergy=turbine_powers_by_direction_thomas2019_new[1,:], FLOWFarm=turbine_powers_by_direction_ff[1,:], Diff=difference_turbines_ffvspe[1,:])
    println(describe(dfffvspe))

    # plot errors
end

# function plantenergy_vs_flowfarm(;nturbines=2)
#     # load plant energy data 
#     turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019 = get_data(journal=false)
#     turbine_powers_by_direction_sowfa, turbine_powers_by_direction_thomas2019_new = get_data(journal=true)

#     # load and run FLOWFarm setup 
#     turbine_powers_by_direction_ff = run_flow_farm(use_local_ti=true, windrose="none", nsamplepoints=100, alpha=1.0, verbose=true, shearfirst=false, filename="../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")

#     println(turbine_powers_by_direction_ff)
#     df = DataFrame(FLOWFarm = turbine_powers_by_direction_ff[:])
#     CSV.write("bp_turb_powers.txt", df)
#     # plot errors
# end