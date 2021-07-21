using Arrow 
using DataFrames
using CSV
using Statistics

function read_les_data(case)
    if case == "low-ti"
        path = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/low_TI/c_000_38turb_D126_lowTI_wd$(10)/"
    else
        path = "/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/data/high_TI/c_000_38turb_D126_hiTI_wd$(10)/"
    end
        filename = "lite_data.ftr"
    table = Arrow.Table(path*filename)
    return DataFrame(table)
end

function get_inflow_u_ave(df::DataFrame, case="high-ti")
    df_plane = df[(df[:, :x].==5.0),:]
    h = unique(df_plane[(df_plane[:, :z].<300.0), :z])
    u = zeros(0)
    for val in h
        tmp_avg = mean(df_plane[(df_plane[:, :z].==val),:u])
        println(val, tmp_avg)
        push!(u, tmp_avg)
    end
    
    df_out = DataFrame(h=h, u=u)

    CSV.write("wind-shear-les-$case.txt", df_out)
end

function parse_ftr(case)
    df = read_les_data(case)
    get_inflow_u_ave(df, case)
end

