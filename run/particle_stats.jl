
include("../src/build.jl")

using .Threads

function main(ARGS)
    println(" ")
    println("=========")
    println(ARGS)
    println("=========")
    println(" ")
    @time compute(ARGS)
end

function compute(ARGS::Vector{String})
    case = ARGS[1]; edir = ARGS[2]
    ti   = parse(Int, ARGS[3]); tf = parse(Int, ARGS[4])
    n = parse(Int, ARGS[5])

    # ps = sample(get_parts(edir*"/particles/", ti), 25_000, replace=false)
    # @time uu, s = uu_cond_r(ps, 128, (1, false))
    # uu = [mean([uu[1,1,i], uu[2,2,i], uu[3,3,i]]) for i in 1:128]
    # s = mean(diag(s))
    # writedlm("outs/uu/bf_"*case*".txt", uu)
    # writedlm("outs/uu/bf_"*case*"_s.txt", s)

    #  ps = sample(get_parts(edir*"/particles/", ti), 25_000, replace=false)
    #  @time uu, s = vv_cond_r(ps, 128, (1, false))
    #  uu = [mean([uu[1,1,i], uu[2,2,i], uu[3,3,i]]) for i in 1:128]
    #  s = mean(diag(s))
    #  writedlm("outs/vv/bf_"*case*".txt", uu)
    #  writedlm("outs/vv/bf_"*case*"_s.txt", s)

    #  println("Getting RDF...")
    #  @time prdf(case, edir, ti, tf, n)
    #  println("RDF done")
    #  println(" ")

    # println("Getting R_uu...")
    # @time puu(case, edir, ti, tf, n)
    # println("R_uu done")
    # println(" ")

    println("Getting R_vv...")
    @time pvv(case, edir, ti, tf, n)
    println("R_vv done")
    println(" ")
    return 
end

function prdf(case::String, edir::String, ti::Int64, tf::Int64, n::Int64)
    rdf = zeros(Float64, (nthreads(), 128)); c = zeros(Float64, nthreads()); drv = zeros(Float64, nthreads())
    @threads for step in ti:tf
        ps = get_parts(edir*"/particles/", step)
        hit= get_grid(edir, step, n, 6.2832)
        dr, rdft = pic_rdf(sample(ps, 25_000, replace=false), hit, 0.6, 128; ncells=16)
        rdf[threadid(), :] += rdft
        c[threadid()] += 1.
        drv[threadid()] += dr
    end
    rdf = reshape(sum(rdf, dims=1), :)/sum(c)
    dr = sum(drv)/sum(c)
    writedlm("./outs/rdf/"*case*".txt", rdf)
    writedlm("./outs/rdf/"*case*"_dr.txt", dr)
    return
end

function puu(case::String, edir::String, ti::Int64, tf::Int64, n::Int64)
    uu = zeros(Float64, (nthreads(), 129))
    c = zeros(Float64, nthreads())
    drv = zeros(Float64, nthreads())
    s = zeros(Float64, nthreads())
    @threads for step in ti:tf
        ps = get_parts(edir*"/particles/", step)
        hit= get_grid(edir, step, n, 6.2832)
        dr, uut, st = pic_uu(sample(ps, 1_000, replace=false), hit, 6.2832, 128; ncells=1)
        uu[threadid(), :] += uut
        s[threadid()] += st
        c[threadid()] += 1.
        drv[threadid()] += dr
    end
    uu = reshape(sum(uu, dims=1), :)/sum(c)
    dr = sum(drv)/sum(c)
    s  = sum(s)/sum(c)
    writedlm("./outs/uu/"*case*".txt", uu)
    writedlm("./outs/uu/"*case*"_dr.txt", dr)
    writedlm("./outs/uu/"*case*"_s.txt", s)
    return
end

function pvv(case::String, edir::String, ti::Int64, tf::Int64, n::Int64)
    vv = zeros(Float64, (nthreads(), 128))
    c = zeros(Float64, nthreads())
    drv = zeros(Float64, nthreads())
    s = zeros(Float64, nthreads())
    @threads for step in ti:tf
        ps = get_parts(edir*"/particles/", step)
        hit= get_grid(edir, step, n, 6.2832)
        dr, vvt, st = pic_vv(sample(ps, 15_000, replace=false), hit, 1.5, 128; ncells=16)
        vv[threadid(), :] += vvt
        s[threadid()] += st
        c[threadid()] += 1.
        drv[threadid()] += dr
    end
    vv = reshape(sum(vv, dims=1), :)/sum(c)
    dr = sum(drv)/sum(c)
    s  = sum(s)/sum(c)
    writedlm("./outs/vv/"*case*".txt", vv)
    writedlm("./outs/vv/"*case*"_dr.txt", dr)
    writedlm("./outs/vv/"*case*"_s.txt", s)
    return
end


main(ARGS)
