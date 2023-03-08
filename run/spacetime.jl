include("../src/build.jl")
using .Threads

function parallel_tptt()
    ns=76; nb=128
#    Rp = zeros(Float32, (ns, nb, nthreads()))
#    S = zeros(nthreads())
#    @threads for _ in 1:16
#        R, s = tptt()
#        Rp[:,:,threadid()] .+= R
#        S[threadid()]       += s
#    end
#    Rp = reshape(sum(Rp, dims=3), (ns, nb))./16
#    S  = reshape(sum(S, dims=1), :)./16
    R, S = tptt()
    writedlm("./R.txt", R)
    writedlm("./S.txt", S)
    return
end

function tptt()
    n = 20_000; nsamples=10_000; nb=128; ti=25; tf=100; ns=76
    dir = "/home/maxzog/NGA2s/nga2restart/examples/hit/ensight/HIT/"
    ps  = sample(rand_ps(n), nsamples, ordered=true, replace=false)
    hit = get_grid(dir, ti, 256, 6.2832)    
    update_fld!(ps, hit)

    R, Rc, dr = init_array(hit, nb, ns)

    self=0; sc = 0;

    for t in ti:tf  
	@printf("\nStep : %i", t)
        hit = get_grid(dir, t, 256, 6.2832)
        for i in 1:nsamples
            for j in i:nsamples
                p = ps[i]; q = ps[j]

                if p.id == q.id && t==ti
                    self += dot(p.fld, interp_fld(q.pos, hit))
                    sc   += 1
                end
                r = get_minr(p.pos, q.pos, hit.L)
                rind = floor(Int, r/dr) + 1
                R[t-ti+1 , rind] += dot(p.fld, interp_fld(q.pos, hit))
                Rc[t-ti+1, rind] += 1
            end
        end
    end

    for i in 1:ns
        for j in 1:nb
            if Rc[i,j] != 0
                R[i,j] /= Rc[i,j]
            end
        end
    end
    return R, self/sc
end

function init_array(field::grid, nb::Int64, ns::Int64)
    R = zeros(Float32, (ns, nb)); Rc = zeros(Float32, (ns, nb))
    rmax = norm([field.L/2, field.L/2, field.L/2])
    dr = rmax/nb
    return R, Rc, dr
end

function div_by_count!(arr, c)
    for i in 1:size(arr)[1]
        for j in size(arr)[2]
            if c[i,j] != 0
                arr[i,j] /= c[i,j]
            end
        end
    end
end

function update_fld!(ps::Vector{part}, hit::grid)
    for p in ps
        p.fld[:] = interp_fld(p.pos, hit)
    end
end

function interp_fld(X::Vector{Float32}, field::grid)::Vector{Float32}
    xm = collect(field.Δ/2:field.Δ:field.L-field.Δ/2)
    ym = collect(field.Δ/2:field.Δ:field.L-field.Δ/2)
    zm = collect(field.Δ/2:field.Δ:field.L-field.Δ/2)
    return fluid2part(X, xm, ym, zm, field.U, field.V, field.W)
end

function rand_ps(n::Int64; L=Float32(6.2832))::Vector{part}
    ps = init_parts(n)
    for p in ps
        p.pos[:] = rand(Float32, 3).*L
    end
    return ps
end

parallel_tptt()


