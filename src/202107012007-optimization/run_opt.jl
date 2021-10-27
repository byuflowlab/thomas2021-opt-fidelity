include("202107012007-optimization.jl")

run_optimization_series(1, "low-ti", "sowfa-nrel", "./low-ti-aec-no-wec/", 
"angle-each-circle", wec=true, lspacing=5.0, firstrun=1, verbose=true, plotresults=false)

#run_optimization_series(400, "low-ti", "sowfa-nrel", "./low-ti-aec-no-wec/", 
#"angle-each-circle", wec=true, lspacing=5.0, firstrun=1, verbose=true, plotresults=false)
