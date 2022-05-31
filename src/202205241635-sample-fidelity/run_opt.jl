using ArgParse

include("sample-fidelity-optimization.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--firstrun"
            help = "which layout to start with"
            arg_type = Int
            default = 1
        "--nruns"
            help = "how many layouts to run in this set"
            arg_type = Int
            default = 1
        "--case"
            help = "either high-ti or low-ti"
            default = "low-ti"
            arg_type = String
        "--out-dir"
            help = "where to put generated files"
            arg_type = String
            default = "./"
        "--dir-bins"
            help = "how many directions to use"
            arg_type = Int
            default = 36
        "--nrotorpoints"
        help = "how many rotor sample points to use"
        arg_type = Int
        default = 36
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    run_optimization_series(parsed_args["nruns"], parsed_args["case"], "sowfa-nrel", parsed_args["nrotorpoints"]; outdir=parsed_args["out-dir"], ndirectionbins=parsed_args["dir-bins"], layoutgen="angle-each-circle", wec=true, lspacing=5.0, firstrun=parsed_args["firstrun"], verbose=true, plotresults=false, savehistory=true)

end

main()
