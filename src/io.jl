

abstract type Particle{Float32<:Real} end
abstract type grid{Float32<:Real} end

struct dnsgrid<:grid{Float32}
    n :: Int64
    L :: Float32
    Δ :: Float32
    U :: Array{Float32}
    V :: Array{Float32}
    W :: Array{Float32}
    P :: Array{Float32}
end

struct lesgrid<:grid{Float32}
    n :: Int64
    L :: Float32
    Δ :: Float32
    U :: Array{Float32}
    V :: Array{Float32}
    W :: Array{Float32}
    P :: Array{Float32}
    ν :: Array{Float32}
    λ :: Array{Float32}
    ϵ :: Array{Float32}
end

struct mon
    Reλ    :: Vector{Float32}
    TKE    :: Vector{Float32}
    sgsTKE :: Vector{Float32}
    η      :: Vector{Float32}
    ϵ      :: Vector{Float32}
    urms   :: Vector{Float32}
end

mutable struct part{Float32}<:Particle{Float32}
    id  :: Int64
    pos :: Vector{Float32}
    vel :: Vector{Float32}
    fld :: Vector{Float32}
    uf  :: Vector{Float32}
end

mutable struct part_dns{Float32}<:Particle{Float32}
    id  :: Int64
    pos :: Vector{Float32}
    vel :: Vector{Float32}
    fld :: Vector{Float32}
end

mutable struct ou_part{Float32}<:Particle{Float32}
    id  :: Int64
    pos :: Vector{Float32}
    vel :: Vector{Float32}
    fld :: Vector{Float32}
    a   :: Float32
    b   :: Float32
end

function get_grid(dir::String, step::String, n::Int64, L::Float64)::grid
    """
    Reads NGA2 ensight Eulerian data and returns a grid type containing the velocity and pressure fields
    along with some other relevant parameters (grid size, domain length, number of cells)
    Really this is only for HIT or other cubic domains as the read_vel and read_P funcs assume this shape for arrays
    """
    fnv = dir*"/velocity/velocity." * "0"^(6-length(step)) * step
    fnp = dir*"/pressure/pressure." * "0"^(6-length(step)) * step
    U, V, W = read_vel(fnv, n)
    Δ = Float32(L/n)
    P = read_P(fnp, n)
    return dnsgrid(n, Float32(L), Δ, U, V, W, P)
end

function get_lesgrid(dir::String, step::Int64, n::Int64, L::Float64)::lesgrid
    """
    Reads NGA2 ensight Eulerian data and returns a grid type containing the velocity and pressure fields
    along with some other relevant parameters (grid size, domain length, number of cells)
    Really this is only for HIT or other cubic domains as the read_vel and read_P funcs assume this shape for arrays
    """
    step = string(step)
    fnv = dir*"/velocity/velocity." * "0"^(6-length(step)) * step
    fnp = dir*"/pressure/pressure." * "0"^(6-length(step)) * step
    fnν = dir*"/visc/visc." * "0"^(6-length(step)) * step
    fnλ = dir*"/lambda/lambda." * "0"^(6-length(step)) * step
    fnϵ = dir*"/eps/eps." * "0"^(6-length(step)) * step
    U, V, W = read_vel(fnv, n)
    Δ = Float32(L/n)
    P = read_P(fnp, n)
    ν = read_P(fnν, n)
    λ = read_P(fnλ, n)
    ϵ = read_P(fnϵ, n)
    return lesgrid(n, Float32(L), Δ, U, V, W, P, ν, λ, ϵ)
end

function get_grid(dir::String, step::Int64, n::Int64, L::Float64)::grid
    """
    Reads NGA2 ensight Eulerian data and returns a grid type containing the velocity and pressure fields
    along with some other relevant parameters (grid size, domain length, number of cells)
    Really this is only for HIT or other cubic domains as the read_vel and read_P funcs assume this shape for arrays
    """
    step = string(step)
    fnv = dir*"/velocity/velocity." * "0"^(6-length(step)) * step
    fnp = dir*"/pressure/pressure." * "0"^(6-length(step)) * step
    U, V, W = read_vel(fnv, n)
    Δ = Float32(L/n)
    P = read_P(fnp, n)
    return dnsgrid(n, Float32(L), Δ, U, V, W, P)
end

function get_mon(dir::String)::mon
    rel    = read_mon_col(dir*"/hit", "Re_lambda")
    TKE    = read_mon_col(dir*"/hit", "TKE")
    sgsTKE = read_mon_col(dir*"/hit", "sgsTKE")
    eta    = read_mon_col(dir*"/hit", "eta")
    eps    = read_mon_col(dir*"/hit", "Epsilon")
    rms    = read_mon_col(dir*"/hit", "URMS")
    return mon(rel, TKE, sgsTKE, eta, eps, rms)
end

function get_parts_ou(dir::String, step::Int64)::Vector{ou_part}
    """
    Returns vector of particles with id, position, particle velocity, and fldvel seen
    """
    suf = "."*"0"^(6-length(string(step)))*string(step)
    np = get_npart(dir*"particle"*suf)
    X  = read_pos(dir*"particle"*suf, np)
    U  = read_arr(dir*"fld"*suf, np)
    V  = read_arr(dir*"vel"*suf, np)
    ids= read_vec(dir*"id"*suf, np)
    as = read_vec(dir*"a_crw"*suf, np)
    bs = read_vec(dir*"b_crw"*suf, np)
    ps = Vector{ou_part}(undef, np)
    for i in 1:np
        ps[i] = ou_part(trunc(Int64, ids[i]), X[:,i], V[:,i], U[:,i], as[i], bs[i])
    end
    perm = sortperm([p.id for p in ps])
    return ps[perm]
end


function get_parts(dir::String, step::Int64)
    """
    Returns vector of particles with id, position, particle velocity, and fldvel seen
    """
    suf = "."*"0"^(6-length(string(step)))*string(step)
    np = get_npart(dir*"particle"*suf)
    X  = read_pos(dir*"particle"*suf, np)
    U  = read_arr(dir*"fld"*suf, np)
    V  = read_arr(dir*"vel"*suf, np)
    ids= read_vec(dir*"id"*suf, np)
    Uf = read_arr(dir*"uf"*suf, np)
    ps = Vector{part}(undef, np)
    for i in 1:np
        ps[i] = part(trunc(Int64, ids[i]), X[:,i], V[:,i], U[:,i], Uf[:,i])
    end
    perm = sortperm([p.id for p in ps])
    return ps[perm]
end

function get_parts_dns(dir::String, step::Int64)
    """
    Returns vector of particles with id, position, particle velocity, and fldvel seen
    """
    suf = "."*"0"^(6-length(string(step)))*string(step)
    np = get_npart(dir*"particle"*suf)
    X  = read_pos(dir*"particle"*suf, np)
    U  = read_arr(dir*"fld"*suf, np)
    V  = read_arr(dir*"vel"*suf, np)
    ids= read_vec(dir*"id"*suf, np)
    ps = Vector{part}(undef, np)
    for i in 1:np
        ps[i] = part(trunc(Int64, ids[i]), X[:,i], V[:,i], U[:,i], Float32.([0., 0., 0.]))
    end
    perm = sortperm([p.id for p in ps])
    return ps[perm]
end

function init_parts(n::Int64)::Vector{part}
    """
    init particle vector for n particles
        pos, vel, fld all equal to zero
    """
    ps = [part(i, Float32.(zeros(3)), Float32.(zeros(3)), Float32.(zeros(3)), Float32.(zeros(3))) for i in 1:n]
    return ps
end

function read_P(fn::String, n::Int64)::Array{Float32}
    """
    Reads pressure... what else do you want from me?
    """
    arr = Array{Float32}(undef, (n,n,n))
    io = open(fn)
    skip(io, 244)
    read!(io, arr); close(io)
    return arr
end

function read_vec(fn::String, n::Int32)::Vector{Float32}
    """
    Reads ensight particle data files and outputs a vector of particle scalar data
    The header size can vary (as far as I can tell) so read the particle.****** file first using get_npart()
    to determine how much to skip in the data file
    """
    vec = Vector{Float32}(undef, n)
    io = open(fn)
    skip(io, 80)
    read!(io, vec); close(io)
    return vec
end

function read_arr(fn::String, n::Int32)::Array{Float32}
    """
    Reads ensight particle data files and outputs a (3, npart) array of the vector data
    The header size can vary (as far as I can tell) so read the particle.****** file first using get_npart()
    to determine how much to skip in the data file
    """
    arr = Array{Float32}(undef, (3, n))
    io = open(fn)
    skip(io, 80)
    read!(io, arr); close(io)
    return arr
end

function read_pos(fn::String, n::Int32)::Array{Float32}
    """
    Reads ensight particle geometry files
    outputs (3, npart) array of position data
    """
    arr = Array{Float32}(undef, (3, n))
    io = open(fn, "r")
    skip(io, 244+4*n)
    read!(io, arr)
    close(io); return arr
end

function read_vel(fn::String, n::Int64)::Tuple{Array{Float32}, Array{Float32}, Array{Float32}}
    """
    Reads velocity field
    Returns three arrays (U, V, W)
    """
    U = Array{Float32}(undef, (n,n,n))
    V = Array{Float32}(undef, (n,n,n))
    W = Array{Float32}(undef, (n,n,n))
    io = open(fn)
    skip(io, 244)
    read!(io, U)
    read!(io, V)
    read!(io, W)
    return U, V, W
end

function get_npart(fn::String)::Int32
    """
    Reads the ensight particle.****** files to get the number of particles
    There's probably a more elegant way to do this but hey this works 
    """
    if split(split(fn, '.')[end-1], '/')[end] != "particle"
        print(fn)
        error("[get_npart] Not the correct file for reading npart")
        return 
    end
    io = open(fn)
    skip(io, 240)
    n = Vector{Int32}(undef, 1)
    read!(io, n); close(io)
    return n[1]
end

function man_split(fn::String)
    """
    For reading in NGA2 monitor files
    Returns vector of column headers
    Handles column headers of one or two lines
    """
    io = open(fn, "r")
    ls = readlines(io)
    ncols = lastindex(split(ls[1]))
    cols = Vector{String}(undef, ncols)
    println(length(ls[2]))
    for i in 1:ncols
        if i == 1
            cols[i] = split(ls[1])[i]*" "*strip(ls[2][1:12])
        else
            j = 13 + (i-2)*14
            cols[i] = strip((split(ls[1])[i]*" "*ls[2][j:j+13]))
        end
    end
    return cols, ls
end

function final_step(ls::Vector{String})::Int64
    return parse(Int64, split(ls[end])[1])
end

function read_mon_col(fn::String, col::String)::Vector{Float32}
    """
    Reads NGA2 monitor files
    Outputs vector of [col name] vs. time
    """
    # cols, ls = man_split(fn)
    # println(cols)
    io = open(fn, "r")
    ls = readlines(io)
    cols = strip.(split(ls[1]))
    ind = findall(el->el==col, cols)
    arr = zeros(lastindex(ls))
    for i in lastindex(ls):-1:1
        vals = split(ls[i])
        if lastindex(vals) != lastindex(cols)
            return arr[lastindex(ls)-final_step(ls)+1:end]
        else
            try
                arr[i] = parse(Float32, vals[ind][1])
            catch ArgumentError
                return arr
            end
        end
    end
    return arr[lastindex(ls)-final_step(ls):end]
 end

