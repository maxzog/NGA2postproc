
function fld_autocorr(dir::String, mid::Int64)
    nsteps = [parse(Int, split(f, '.')[end]) for f in glob("id*", dir)]
    nsteps = lastindex(unique(nsteps))
    @assert mid < nsteps
    ac = zeros(Float64, nsteps-mid+1)
    pso = get_parts(dir, nsteps)
    npart = lastindex(pso)
    μf = mean([p.fld for p in pso])
    μs = mean([p.uf for p in pso])
    self = mean([dot(p.fld+p.uf-μs-μf, p.fld+p.uf-μs-μf) for p in pso])
    for i in mid:nsteps
        ps = get_parts(dir, i)
        for j in 1:lastindex(ps)
            ac[i-mid+1] += dot(ps[j].fld+ps[j].uf-μs-μf, pso[j].fld+pso[j].uf-μs-μf)
        end
    end
    return ac/npart#/self
end

function fldr_autocorr(dir::String, mid::Int64)
    nsteps = [parse(Int, split(f, '.')[end]) for f in glob("id*", dir)]
    nsteps = lastindex(unique(nsteps))
    @assert mid < nsteps
    ac = zeros(Float64, nsteps-mid+1)
    pso = get_parts(dir, nsteps)
    npart = lastindex(pso)
    μf = mean([p.fld for p in pso])
    self = mean([dot(p.fld-μf, p.fld-μf) for p in pso])
    for i in mid:nsteps
        ps = get_parts(dir, i)
        for j in 1:lastindex(ps)
            ac[i-mid+1] += dot(ps[j].fld-μf, pso[j].fld-μf)
        end
    end
    return ac/npart#/self
end

function us_autocorr(dir::String, mid::Int64)
    nsteps = [parse(Int, split(f, '.')[end]) for f in glob("id*", dir)]
    nsteps = lastindex(unique(nsteps))
    @assert mid < nsteps
    ac = zeros(Float64, nsteps-mid+1)
    pso = get_parts(dir, nsteps)
    npart = lastindex(pso)
    μ = mean([p.uf for p in pso])
    self = mean([dot(p.uf-μ, p.uf-μ) for p in pso])
    for i in mid:nsteps
        ps = get_parts(dir, i)
        for j in 1:lastindex(ps)
            ac[i-mid+1] += dot(ps[j].uf-μ, pso[j].uf-μ)
        end
    end
    return ac/npart/self
end

function vel_autocorr(dir::String, mid::Int64)
    nsteps = [parse(Int, split(f, '.')[end]) for f in glob("id*", dir)]
    nsteps = lastindex(unique(nsteps))
    @assert mid < nsteps
    ac = zeros(Float64, nsteps-mid+1)
    pso = get_parts(dir, nsteps)
    npart = lastindex(pso)
    μ = mean([p.vel for p in pso])
    self = mean([dot(p.vel-μ, p.vel-μ) for p in pso])
    for i in mid:nsteps
        ps = get_parts(dir, i)
        for j in 1:lastindex(ps)
            ac[i-mid+1] += dot(ps[j].vel-μ, pso[j].vel-μ)
        end
    end
    return ac/npart#/self
end

function uu_twotime()
    
end
