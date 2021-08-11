using FLOWFarm; const ff=FLOWFarm 
using DataFrames
using Statistics
include("../202106291409-les-data-parsing/202106291410-parse-les-ftr.jl")
#read_les_data, get_inflow_u_ave, parse_ftr

function get_uprime_at_h(df, h)

    # subset of flowfield for free-stream points at hub height 
    data = df[(df.z .== h) .& (df.x .< 500.0), :]
    
    # actual velocity 
    u1d = data.u #hypot.(data.u, data.v, data.w)

    # mean velocity
    ubar = mean(u1d)

    # turbulent fluctuations 
    uprime = u1d .- ubar

    return uprime, ubar
end

function get_average_ti(case="high-ti")
    # based on https://www.mit.edu/course/1/1.061/www/dream/SEVEN/SEVENTHEORY.PDF

    # load dataframe of flow field
    df = read_les_data(case)

    uprime = []
    ubar = []
    urms = []
    # for h = 25:10:155
    for h = 75:10:105
        uprimet, ubart = get_uprime_at_h(df, h)
        push!(uprime, uprimet)
        push!(ubar, ubart)

        # turbulence strength
        # urmst = sqrt((1.0/length(uprimet))*sum(uprimet.^2))
        # urmst = sqrt(sum(abs.(uprimet))/length(uprimet))/ubart
        urmst = mean(sqrt.(uprimet.^2)/ubart)
        push!(urms, urmst)
    end
    println(urms)
    println("\nMean TI = $(round.(mean(urms), digits=3))\nStd. TI = $(round.(std(urms), digits=3))\nWS = $(round.(mean(ubar), digits=3))")

end