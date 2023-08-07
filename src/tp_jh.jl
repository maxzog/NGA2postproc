
function pic_uu_lt(ps::Vector{part_dns}, L::Float32, rmax::Float64, nb::Int64; ncells=16)
    re_id!(ps)
    Δ = L / ncells
    npic, ipic = part_grid(ps, Δ, ncells)
    np = lastindex(ps)
    @assert sum(npic) == np
    #println("Warning: Overwriting all p.uf to zero")
    #for p in ps
    #    p.uf = Float32.([0., 0., 0.])
    #end

    no = floor(Int64, rmax/Δ)
    if no == 0
        rmax = norm([Δ, Δ, Δ])
    else
        rmax = (2*no+1)*norm([Δ, Δ, Δ])
    end
    dr = rmax/nb
    uul = zeros(Float64, nb+1); rv = 0:dr:rmax
    uut = zeros(Float64, nb+1); rv = 0:dr:rmax
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
                            r = get_minr(p.pos, q.pos, L)
                            ir = floor(Int, r/dr) + 1
                            if p.id != q.id
                                rl, rt = par_perp_u(p, q)
                                uul[ir] += dot(p.fld,rl)*dot(q.fld,rl)
                                uut[ir] += dot(p.fld,rt)*dot(q.fld,rt) 
                                c[ir]  += 1
                            end
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
    uul = uul./(3*c); uut = uut./(3*c); s = s/(3*sc)
    return dr, uul, uut, s
end

function pic_struct_lt(ps::Vector{part_dns}, L::Float32, rmax::Float64, nb::Int64; ncells=16)
    re_id!(ps)
    Δ = L / ncells
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
    uul = zeros(Float64, nb+1); rv = 0:dr:rmax
    uut = zeros(Float64, nb+1); rv = 0:dr:rmax
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
                            r = get_minr(p.pos, q.pos, L)
                            ir = floor(Int, r/dr) + 1
                            if p.id != q.id
                                rl, rt = par_perp_u(p, q)
                                url = dot(q.fld,rl)-dot(p.fld,rl)
                                urt = dot(q.fld,rt)-dot(p.fld,rt)
                                uul[ir] += url*url
                                uut[ir] += urt*urt
                                c[ir]  += 1
                            end
                        end
                        if ni != 0 && ps[ni].id == p.id
                            s += dot(ps[ni].fld-p.fld,ps[ni].fld-p.fld)
                            sc += 1
                        end
                    end
                end
            end
        end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    uul = uul./c; uut = uut./c; s = s/(3*sc)
    return dr, uul, uut, s
end

function par_perp_gen(r::Vector{Float32}, L::Float32, vel::Vector{Float32})
    for i in 1:3
       if r[i] > L/2
          r[i] -= L
       end
       if r[i] < -L/2
          r[i] += L
       end
    end
    rll = r/norm(r)
    rt1 = cross(rll,vel) / norm(cross(rll,vel))
    rt2 = cross(rt1,rll) / norm(cross(rt1,rll))
    @assert norm(rll)≈1.
    @assert norm(rt1)≈1.
    @assert norm(rt2)≈1.
    return rll, rt2
end

function par_perp_u(p::part, q::part, L::Float32)
    r = q.pos-p.pos
    for i in 1:3
       if r[i] > L/2
          r[i] -= L
       end
       if r[i] < -L/2
          r[i] += L
       end
    end
    rll = r/norm(r)
    rt1 = cross(rll,p.fld+p.uf) / norm(cross(rll,p.fld+p.uf))
    rt2 = cross(rt1,rll) / norm(cross(rt1,rll))
    @assert norm(rll)≈1.
    @assert norm(rt1)≈1.
    @assert norm(rt2)≈1.
    return rll, rt2
end

function par_perp_u(p::part_dns, q::part_dns, L::Float32)
   r = q.pos-p.pos
   for i in 1:3
      if r[i] > L/2
         r[i] -= L
      end
      if r[i] < -L/2
         r[i] += L
      end
   end
   rll = r/norm(r)
   rt1 = cross(rll,p.fld) / norm(cross(rll,p.fld))
   rt2 = cross(rt1,rll) / norm(cross(rt1,rll))
   @assert norm(rll)≈1.
   @assert norm(rt1)≈1.
   @assert norm(rt2)≈1.
   return rll, rt2
end


function pic_sff(ps::Vector{part}, field::grid, rmax::Float64, nb::Int64; ncells=16)
    re_id!(ps)
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
                            uu[ir] += dot(p.fld-q.fld, p.fld-q.fld)
                            c[ir]  += 1
                        end
                        if ni != 0 && ps[ni].id == p.id
                            s += dot(p.fld-ps[ni].fld, p.fld-ps[ni].fld)
                            sc += 1
                        end
                    end
                end
            end
        end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    uu = uu./(3*c); s = s/(3*sc)
    return dr, uu, s
end


function pic_uu(ps::Vector{part_dns}, field::grid, rmax::Float64, nb::Int64; ncells=16)
    re_id!(ps)
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
    return dr, uu, s
end
function re_id!(ps::Vector{part_dns})
    for i in 1:lastindex(ps)
        ps[i].id = i
    end
end

function re_id!(ps::Vector{part})
    for i in 1:lastindex(ps)
        ps[i].id = i
    end
end

function part_grid(ps::Vector{part_dns}, Δ::Float32, ncells::Int64)::Tuple{Array{Int64}, Array{Int64}}
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

function get_ijk(p::part_dns, n::Int64, Δ::Float32)::Tuple{Int64, Int64, Int64}
    i, j, k = trunc.(Int64, ceil.(p.pos ./ Δ))
    @assert i <= n; @assert i>=1
    @assert j <= n; @assert j>=1
    @assert k <= n; @assert k>=1
    return (i,j,k)
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

function uu_cond_r(ps::Vector{part_dns}, nbins::Int64, L::Float32)
    uu = zeros(Float32, nbins+1)
    rmax = norm([L/2, L/2, L/2])
    dr = rmax / nbins
    counts = zeros(Int64, nbins+1)
    scounts = 0
    muf = mean([p.fld for p in ps])
    self = 0.0
    for i ∈ 1:lastindex(ps), j ∈ i:lastindex(ps)
        r = norm(ps[i].pos-ps[j].pos)
        ir = floor(Int, r/dr) + 1
        ir > lastindex(uu) ? check = false : check = true
        if check
           uu[ir] += dot(ps[i].fld-muf, ps[j].fld-muf)
           counts[ir] += 1
        end
        if i == j
            self += dot(ps[i].fld-muf, ps[j].fld-muf)
            scounts += 1
        end
    end
    self /= scounts
    for i in 1:lastindex(counts)
        if counts[i] != 0
            uu[i] = uu[i] / counts[i]
        end
    end
    return dr, uu, self
end

function uu_lt_cond_r(ps::Vector{part_dns}, nb::Int64, L::Float32)
    re_id!(ps)
    rmax = norm([L/2 L/2 L/2])
    dr = rmax/nb
    rv = 0:dr:rmax
    uul = zeros(Float64, nb+1)
    uut = zeros(Float64, nb+1)
    c  = zeros(Int, nb+1); sc = 0; s = 0.

    for i in 1:lastindex(ps), j in i:lastindex(ps)
       p = ps[i]; q = ps[j]
       # r = get_minr(p.pos, q.pos, L)
       r = norm(q.pos-p.pos)
       ir = floor(Int, r/dr) + 1
       ir > lastindex(uul) ? check = false : check = true
       if p.id != q.id && check
           rl, rt = par_perp_u(p, q)
           uul[ir] += dot(p.fld,rl)*dot(q.fld,rl)
           uut[ir] += dot(p.fld,rt)*dot(q.fld,rt) 
           c[ir]  += 1
       end
       if q.id == p.id
           s += dot(p.fld, q.fld)
           sc += 1
       end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    uul = uul./(3*c); uut = uut./(3*c); s = s/(3*sc)
    return dr, uul, uut, s
end

function sf_lt_cond_r(ps::Vector{part_dns}, nb::Int64, L::Float32)
    re_id!(ps)
    rmax = norm([L/2 L/2 L/2])
    dr = rmax/nb
    rv = 0:dr:rmax
    sfl = zeros(Float64, nb+1)
    sft = zeros(Float64, nb+1) 
    c  = zeros(Int, nb+1); sc = 0; s = 0.
    for i in 1:lastindex(ps)
       for j in i:lastindex(ps)
         p = ps[i]; q = ps[j]
         # r = get_minr(p.pos, q.pos, L)
         r = norm(q.pos-p.pos)
         ir = floor(Int, r/dr) + 1
         ir > lastindex(sfl) ? check = false : check = true
         if p.id != q.id && check
             rl, rt = par_perp_u(p, q, L)
             sfl[ir] += dot(q.fld - p.fld, rl) * dot(q.fld - p.fld,rl)
             sft[ir] += dot(q.fld - p.fld, rt) * dot(q.fld - p.fld,rt) 
             c[ir]  += 1
         end
         if q.id == p.id
             s += dot(q.fld - p.fld, q.fld - p.fld)
             sc += 1
         end
      end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    sfl = sfl./c; sft = sft./c; s = s/(3*sc)
    return dr, sfl, sft, s
end

function general_structure(pos::Vector{Vector{Float32}}, vel::Vector{Vector{Float32}}, nb::Int64, L::Float32)
    rmax = norm([L/2 L/2 L/2])
    dr = rmax/nb
    rv = 0:dr:rmax
    sfl = zeros(Float64, nb+1)
    sft = zeros(Float64, nb+1) 
    c  = zeros(Int, nb+1); sc = 0; s = 0.
    for i in 1:lastindex(pos)
       for j in i:lastindex(pos)
         # r = get_minr(p.pos, q.pos, L)
         r = norm(pos[j]-pos[i])
         ir = floor(Int, r/dr) + 1
         ir > lastindex(sfl) ? check = false : check = true
         if i != j && check
            rl, rt = par_perp_gen(pos[j]-pos[i], L, vel[i])
            sfl[ir] += dot(vel[j]-vel[i], rl) * dot(vel[j]-vel[i],rl)
            sft[ir] += dot(vel[j]-vel[i], rt) * dot(vel[j]-vel[i],rt) 
            c[ir]  += 1
         end
         if i == j
            s += dot(vel[j]-vel[i], vel[j]-vel[i])
            sc += 1
         end
      end
    end
    c = [c[i]==0 ? c[i]=1 : c[i]=c[i] for i in 1:lastindex(c)]
    sfl = sfl./c; sft = sft./c; s = s/(3*sc)
    return dr, sfl, sft, s
end
