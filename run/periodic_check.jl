
include("../src/build.jl")

function main()
    np = 25_000; L = Float32(6.2832)
    ps = rand_ps(np)
    for i in 1:np
        p = ps[i]
        for j in i+1:np
            q = ps[j]
            r1 = get_minr(p.pos, q.pos, L)
            r2 = per_r(p.pos, q.pos, L)
            @assert r1 ≈ r2
            println(r1-r2)
        end
    end
end

function per_r(x1::Vector{Float32}, x2::Vector{Float32}, L::Float32)
    r = zeros(Float32, lastindex(x1))
    for i in 1:lastindex(x1)
        r[i] = x2[i]-x1[i] - 2*((x2[i]-x1[i])%(L/2))
        r[i] < 0 ? r[i] += L/2 : nothing
    end
    return norm(r)
end

function rand_ps(n::Int64; L=Float32(6.2832))::Vector{part}
    ps = init_parts(n)
    for p in ps
        p.pos[:] = rand(Float32, 3).*L
    end
    return ps
end


main()