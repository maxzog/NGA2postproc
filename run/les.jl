#! usr/bin/env julia

include("../src/build.jl")

function main(ARGS)::nothing
    case=ARGS[1]; fn = ARGS[2]; step = parse(Int, ARGS[3])
    ps = sample(get_parts(fn*"/particles/", step), 
                20_000, # number of particles to be used in sample
                replace=false, ordered=true)
    hit= get_lesgrid(fn, step, 32, 6.2832)
    rv, rdf = pic_rdf(ps, hit, 1.0, 100; ncells=32)
    writedlm("./rdfs/"*case*".txt", rdf)
    writedlm("./rdfs/"*case*"_rv.txt", rv)
    rv, rdf = get_rdf(ps, 100, Float32(6.2832), 3, true)
    writedlm("./rdfs/"*case*"old.txt", rdf)
    writedlm("./rdfs/"*case*"_rv_old.txt", rv)
    return nothing
end

main(ARGS)