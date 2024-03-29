
function fftcorrgrid(U::Array{Float32}, comp::Int64)
    sumd = [1,2,3]
    deleteat!(sumd, findall(x->x==comp,sumd))
    uh = fft(U, [comp])
    uhs = conj(uh)
    u = sum(real.(ifft(uh.*uhs, [comp])), dims = sumd) / lastindex(U[:,1,1])^3
    return dropdims(u, dims = Tuple(sumd))[1:trunc(Int64, lastindex(u)/2)]
end

function fftcorrgrid_twotime(U::Array{Float32}, Us::Array{Float32}, comp::Int64)
    sumd = [1,2,3]
    deleteat!(sumd, findall(x->x==comp,sumd))
    uh = fft(U, [comp])
    uhs = conj(fft(Us, [comp]))
    u = sum(real.(ifft(uh.*uhs, [comp])), dims = sumd) / lastindex(U[:,1,1])^3
    return dropdims(u, dims = Tuple(sumd))[1:trunc(Int64, lastindex(u)/2)]
end

function R11(hit::grid)
    corr = (fftcorrgrid(hit.U, 1) + fftcorrgrid(hit.V, 2) + fftcorrgrid(hit.W, 3))./3
    return corr
end

function R22(hit::grid)
    corr = (fftcorrgrid(hit.U, 2) + fftcorrgrid(hit.V, 3) + fftcorrgrid(hit.W, 1))./3
    return corr
end

function uvar(U::Array{Float32})
    up = U .- mean(U)
    return sum(up.*up) ./ lastindex(up[:,1,1])^3
end

function L11(hit::grid)
    R = R11(hit) 
    R ./= R[1]
    return ∫(R, hit.Δ)
end

function L22(hit::grid)
    R = R22(hit) 
    R ./= R[1]
    return ∫(R, hit.Δ)
end

function ∫(f::Vector{Float32}, dx::Float32)
    return sum(f.*dx) 
end

function get_rdf(ps::Vector{part}, nbins::Int64, L::Float32, dim::Int64, nga::Bool)
    npart = length(ps)
    xs = [p.pos for p in ps]
    rmax = norm([L/2 for i in 1:dim])
    dr = rmax / nbins
    rv = 0:dr:rmax
    rdf = zeros(nbins)
    if nga
        N = npart*(npart-1) / 2 
    else
        N = npart - 1
    end
    dim == 2 ? V = L^2 : V = L^3 
    ρ = N/V
    nga ? ll = npart : ll = 1
    for i ∈ 1:ll
        for j ∈ i:npart
            if i != j
                r = get_minr(xs[i], xs[j], L)
                rind = floor(Int, r/dr) + 1
                try
                    rdf[rind] += 1
                catch BoundsError
                    if r ≈ rmax
                        rdf[end] += 1
                    else
                        println("WARNING: Periodicity failure. Breaking.")
                        println("r = ", r)
                        println("Particle ids: ", (i,j))
                        return
                    end
                end
            end
        end
    end
    
    for (i, bin) in enumerate(rdf)
        r = (rv[i] + rv[i+1])/2
        den = 4/3*π*(rv[i+1]^3 - rv[i]^3)

        # The rdf still gets funky past rmax/2... Tried adding corrections from online doc, not working

        # if r < rmax/2
        #     den = 4*π*r^2 * dr
        # elseif r < sqrt(2)*rmax/2
        #     den = 2*π*r*(3-4*r)
        # else
        #     den = 2*r*(3*π -12*f₁(r) + f₂(r))
        # end

        if dim == 2
            den = π * (rv[i+1]^2 - rv[i]^2)
        end
        rdf[i] = bin/(den*ρ)
    end
    return dr/2:dr:rmax-dr/2, rdf
end

function uu_cond_r(psn::Vector{part}, nbins::Int64, subset::Tuple{Int64, Bool};
    L = Float32(1.), Stk = Float32(1), nga = false)
    npart = lastindex(psn)
    if subset[2] && subset[1] != npart
        psn = sample(psn, subset[1], replace = false, ordered = true)
    end
    npart = lastindex(psn)
    usn = [p.fld for p in psn]
    uu = zeros(3,3,nbins)
    rmax = Float32(norm([L/2, L/2, L/2]))
    dr = rmax / nbins
    counts = zeros(nbins)
    scounts = 0
    mu = mean([p.fld for p in psn])
    self = zeros(3,3)
    if nga
        for i ∈ 1:npart, j ∈ i:npart
            r = get_minr(psn[i].pos, psn[j].pos, L)
            rind = floor(Int, r/dr) + 1
            uu[1:3, 1:3, rind] += (usn[i] - mu) * (usn[j] - mu)'
            counts[rind] += 1
            if i == j
                self += (usn[i] - mu) * (usn[j] - mu)'
                scounts += 1
            end
        end
    else
        for i ∈ 1:1, j ∈ 1:npart
            r = get_minr(psn[i].pos, psn[j].pos, L)
            rind = floor(Int, r/dr) + 1
            uu[1:3, 1:3, rind] += (usn[i] - mu) * (usn[j] - mu)'
            counts[rind] += 1
            if i == j
                self += (usn[i] - mu) * (usn[j] - mu)'
                scounts += 1
            end
        end
    end
    self /= scounts
    for i in 1:lastindex(counts)
        if counts[i] != 0
            uu[1:3,1:3,i] = uu[1:3,1:3,i] / counts[i]
        end
    end
    return uu, self
end