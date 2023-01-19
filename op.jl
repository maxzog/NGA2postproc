
function variance(arr::Vector{Vector{Float32}})::Float32
    μ = mean(arr)
    var = mean([dot(x-μ, x-μ)/3 for x in arr])
    return var
end

function autocorrfft(vec::Vector{Float32})::Vector{Float32}
    n = lastindex(vec)
    vec = [zeros(Float32, lastindex(vec)); vec; zeros(Float32, lastindex(vec))]
    uf = fft(vec)
    ufs = conj(uf)
    R = real.(ifft(uf.*ufs))
    return R[1:n]/(3*lastindex(R)/2)
end

function average_monitor(fn::String, col::String, ti::Int64; tf=-1)
    """
    Averages monitor file data in time over [ti, tf] where ti/tf are timestep numbers
    if no tf is supplied, average is taken from ti to the finale timestep
    """
    dat = read_mon_col(fn, col)
    if ti > tf
        tf = lastindex(dat)
    end 
    return mean(dat[ti:tf])
end
