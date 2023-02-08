
function pic_rdf(psin::Vector{part}, field::grid, rmax::Float64, nb::Int64; ncells=16)
    ps = re_id(psin)
    Δ = field.L / ncells
    npic, ipic = part_grid(ps, Δ, ncells)
    np = lastindex(ps); N = np*(np-1)
    @assert sum(npic) == np

    no = floor(Int64, rmax/Δ)
    if no == 0
        rmax = norm([Δ, Δ, Δ])
    else
        rmax = (2*no+1)*norm([Δ, Δ, Δ])
    end
    dr = rmax/nb
    rdf = zeros(Float64, nb); rv = 0:dr:rmax

    V = (field.L)^3
    ρ = N / V
    for p in ps
        ip1, jp1, kp1 = get_ijk(p, ncells, Δ)
        for i in ip1-no:ip1+no
            for j in jp1-no:jp1+no
                for k in kp1-no:kp1+no
                    ii, jj, kk = periodic_inds([i,j,k], ncells)
                    for ni in ipic[:,ii,jj,kk]
                        if ni != 0 && ps[ni].id != p.id
                            q = ps[ni]
                            r = get_minr(p.pos, q.pos, field.L)
                            ir = floor(Int, r/dr) + 1
                            rdf[ir] += 1.
                        end
                    end
                end
            end
        end
    end

    # normalize
    for (i, bin) in enumerate(rdf)
        den = 4/3*π*(rv[i+1]^3 - rv[i]^3)
        rdf[i] = bin/(den*ρ)
    end
    return dr/2:dr:rmax-dr/2, rdf
end

function pic_uu(ps::Vector{part}, field::grid, rmax::Float64, nb::Int64; ncells=16)
    ps = re_id(ps)
    Δ = field.L / ncells
    npic, ipic = part_grid(ps, Δ, ncells)
    np = lastindex(ps)
    @assert sum(npic) == np

    no = floor(Int64, rmax/Δ)
    if no == 0
        rmax = norm([Δ, Δ, Δ])
    else
        rmax = (2*no+1)*norm([Δ, Δ, Δ])
    end
    dr = rmax/nb
    uu = zeros(Float64, nb+1); rv = 0:dr:rmax
    c  = zeros(Int, nb+1); sc = 0; s = 0.

    for p in ps
        ip1, jp1, kp1 = get_ijk(p, ncells, Δ)
        for i in ip1-no:ip1+no
            for j in jp1-no:jp1+no
                for k in kp1-no:kp1+no
                    ii, jj, kk = periodic_inds([i,j,k], ncells)
                    for ni in ipic[:,ii,jj,kk]
                        if ni != 0
                            q = ps[ni]
                            r = get_minr(p.pos, q.pos, field.L)
                            ir = floor(Int, r/dr) + 1
                            uu[ir] += dot(p.fld, q.fld)
                            c[ir]  += 1
                        end
                        if ni != 0 && ps[ni].id == p.id
                            s += dot(p.fld, ps[ni].fld)
                            sc += 1
                        end
                    end
                end
            end
        end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    uu = uu./(3*c); s = s/(3*sc)
    return rv, uu, s
end

function pic_vv(ps::Vector{part}, field::grid, rmax::Float64, nb::Int64; ncells=16)
    ps = re_id(ps)
    Δ = field.L / ncells
    npic, ipic = part_grid(ps, Δ, ncells)
    np = lastindex(ps)
    @assert sum(npic) == np

    no = floor(Int64, rmax/Δ)
    if no == 0
        rmax = norm([Δ, Δ, Δ])
    else
        rmax = (2*no+1)*norm([Δ, Δ, Δ])
    end
    dr = rmax/nb
    vv = zeros(Float64, nb+1); rv = 0:dr:rmax
    c  = zeros(Int, nb+1); sc = 0; s = 0.

    for p in ps
        ip1, jp1, kp1 = get_ijk(p, ncells, Δ)
        for i in ip1-no:ip1+no
            for j in jp1-no:jp1+no
                for k in kp1-no:kp1+no
                    ii, jj, kk = periodic_inds([i,j,k], ncells)
                    for ni in ipic[:,ii,jj,kk]
                        if ni != 0
                            q = ps[ni]
                            r = get_minr(p.pos, q.pos, field.L)
                            ir = floor(Int, r/dr) + 1
                            vv[ir] += dot(p.vel, q.vel)
                            c[ir]  += 1
                        end
                        if ni != 0 && ps[ni].id == p.id
                            s += dot(p.vel, ps[ni].vel)
                            sc += 1
                        end
                    end
                end
            end
        end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    vv = vv./(3*c); s = s/(3*sc)
    return rv, vv, s
end

function re_id(ps::Vector{part})
    for i in 1:lastindex(ps)
        ps[i].id = i
    end
    return ps
end

function part_grid(ps::Vector{part}, Δ::Float32, ncells::Int64)::Tuple{Array{Int64}, Array{Int64}}
    npic = zeros(Int64, (ncells, ncells, ncells))
    for p in ps
        i, j, k = get_ijk(p, ncells, Δ)
        npic[i,j,k] += 1
    end
    maxnpic = maximum(Int64, npic)
    ipic = zeros(Int64, (maxnpic, ncells, ncells, ncells))
    npic .*= 0
    for p in ps
        i, j, k = get_ijk(p, ncells, Δ)
        npic[i,j,k] += 1
        ipic[npic[i,j,k], i, j ,k] = p.id
    end
    return npic, ipic
end

function get_ijk(p::part, n::Int64, Δ::Float32)::Tuple{Int64, Int64, Int64}
    i, j, k = trunc.(Int64, ceil.(p.pos ./ Δ))
    @assert i <= n; @assert i>=1
    @assert j <= n; @assert j>=1
    @assert k <= n; @assert k>=1
    return (i,j,k)
end

function periodic_inds(ijk::Vector{Int64}, n::Int64)
    for d in 1:lastindex(ijk)
        ijk[d]<1 ? ijk[d]+=n : nothing
        ijk[d]>n ? ijk[d]-=n : nothing
    end
    return ijk
end

function correct_norm(r::Float64)::Float64
    if 2*r <= 1
        den = 4*π*r^2
    elseif 2*r < sqrt(2)
        den = 2*π*r*(3-4*r)
    else
        den = 2*r*(3*π -12*f₁(r) + f₂(r))
    end
    return den
end

function f₁(r::Float64)
    return atan(sqrt(4*r^2-1))
end

function f₂(r::Float64)
    return 8*r*atan(2*r*(4*r^2-3)/(sqrt(4*r^2-2)*(4*r^2+1)))
end