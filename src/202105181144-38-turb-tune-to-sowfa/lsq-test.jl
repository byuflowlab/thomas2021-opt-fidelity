using LsqFit 

function obj_func(x, p)
    y = p[1] .+ x.*p[2] 
    ers = (y .- [1., 2., 3., 4., 10., 6., 7., 8., 9., 10.])./50.0
    println("p: $p, error-sum: $((sum(ers.^2))), error-max: $(maximum(abs.(ers)))")
    return y
end

y = [1., 2., 3., 4., 10., 6., 7., 8., 9., 10.]
fitres = curve_fit(obj_func, collect(1.0:10.0), y, [1.0, 1.0])

println(fitres.param)