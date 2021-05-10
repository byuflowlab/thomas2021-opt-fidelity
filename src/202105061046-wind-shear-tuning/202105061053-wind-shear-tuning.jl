using FLOWFarm; const ff=FLOWFarm
using CSV
using DataFrames
using Plots
using StatsPlots
using LsqFit

function obj_func(h, p)

    # constants
    reference_height = 80.0
    reference_speed = 8.0
    ground_height = 0.0

    # get number of points for fitting
    npoints = length(h)

    # initialize shear model 
    shear_model = ff.PowerLawWindShear(p[1])

    # calc v by model
    v_calc = zeros(npoints)
    for i in 1:npoints
        v_calc[i] = ff.adjust_for_wind_shear(h[i], reference_speed, reference_height, ground_height, shear_model)
    end
    return v_calc
end

# set data file name for LES data
datafile = "../inputfiles/wind-shear-les.txt"
# load data to data frame
df = DataFrame(CSV.File(datafile, header=0, datarow=2))
# rename columns
rename!(df,:Column1 => :h,:Column2 => :s)

# optimize fit
initial_shear = 0.1
fit = curve_fit(obj_func, df.h, df.s, [initial_shear])

# final_shear = 0.35
final_shear = fit.param[1]

# set resolution of model data
res = 250
# initialize output
v_calc = zeros(res)
# initialize heights (input)
h_calc = collect(1:res)
# define shear model instance
shear_model = ff.PowerLawWindShear(final_shear)

# calculate speeds at each height
for i in 1:res
    v_calc[i] = ff.adjust_for_wind_shear(h_calc[i], df.s[4], df.h[4], 0.0, shear_model)
end

# put model results in data frame
df2 = DataFrame([h_calc, v_calc], :auto)
rename!(df2,:x1 => :h,:x2 => :s)

# plot LES data and model results
p = scatter(df.s, df.h, label="LES")
plot!(df2.s, df2.h, label="Model")

# print optimized shear value
println("optimized shear: ", final_shear)
println("optimized shear approx.: ", round(final_shear, digits=2))