
function get_minr(X1::Vector{Float32}, X2::Vector{Float32}, perL::Float32; eps = 1e-10)::Float32
    """
    Minimum particle separation distance assuming cubic, periodic domain
    """
    r = zeros(length(X1))
    for i ∈ 1:length(X1)
        dx = X2[i] - X1[i]
        abs(dx) > abs(dx - perL) ? dx -= perL : nothing
        abs(dx) > abs(dx + perL) ? dx += perL : nothing
        r[i] = dx
    end
    return norm(r)
end

function spectral(U, V, W)
    p = plan_fft(U)
    Uk = p*U; Vk = p*V; Wk = p*W
    return Uk, Vk, Wk
end

function ispectral(Uk, Vk, Wk)
    ip = plan_ifft(Uk)
    U = ip*Uk; V = ip*Vk; W = ip*Wk
    return U, V, W
end

function sharpspec(Uk, Vk, Wk, n::Int64, kv::Vector{Int64}, Δ::Float32)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        kk = norm([kv[i] kv[j] kv[k]])
        heavy = π/Δ - kk
        if heavy <= 0
            Uk[i,j,k] = 0
            Vk[i,j,k] = 0
            Wk[i,j,k] = 0
        end
    end
    return Uk, Vk, Wk
end

function spectralfilter(U, V, W, ftype::String, Δ::Float32)
    Uk, Vk, Wk = spectral(U, V, W)
    n = trunc(Int, size(U)[1]/2)
    kv = fftshift(-n:n-1)
    if ftype == "sharpspec"
        Uk, Vk, Wk = sharpspec(Uk, Vk, Wk, 2*n, kv, Δ)
    end
    U, V, W = ispectral(Uk, Vk, Wk)
    return U, V, W
end

function coarsen(U::Array{Float32}, n::Int64)
    @assert n%2 == 0
    nc = trunc(Int, n/2)
    Uc = Array{Float32}(undef, (nc, nc, nc))
    for i in 1:nc
        for j in 1:nc
            for k in 1:nc
                Uc[i,j,k] = mean(U[2*i-1:2*i,2*j-1:2*j,2*k-1:2*k])
            end
        end
    end
    return Uc
end

function box_ind(i::Int64, Δ::Int64, n::Int64)#::Tuple{Int64, Int64}
    is = floor(Int, i-Δ/2); ie = ceil(Int, i+Δ/2)
    ir = collect(is+1:ie-1)
    for j in 1:lastindex(ir)
        ir[j]<1 ? ir[j]+=n : nothing
        ir[j]>n ? ir[j]-=n : nothing
    end
    return ir
end

function box_filt(U::Array{Float32}, Δ::Int64, n::Int64)
    Uf = zeros((n,n,n))
    flen = 0
    for i in 1:n
        ir = box_ind(i,Δ,n)
        for j in 1:n
            jr = box_ind(j,Δ,n)
            for k in 1:n
                kr = box_ind(k,Δ,n)
                cut = U[ir,jr,kr]; flen = cut()
                Δic = 1/size(cut)[1]^3
                @assert sum([Δic for _ in 1:lastindex(cut)]) ≈ 1
                Uf[i,j,k] = Δic * sum(cut)
            end
        end
    end
    return Uf
end

function LES(U, V, W, ftype::String, δ::Int64, L::Float32; eps = 1e-6)
    n = size(U)[1]
    Δ = δ * L / n

    if ftype == "box_conv"
        Ufilt = box_filt(U, δ, n)
        Vfilt = box_filt(V, δ, n)
        Wfilt = box_filt(W, δ, n)
    elseif ftype == "coarsen"
        Ufilt = coarsen(U, n)
        Vfilt = coarsen(V, n)
        Wfilt = coarsen(W, n)
    elseif ftype == "sharpspec"
        Ufilt, Vfilt, Wfilt = spectralfilter(U, V, W, ftype, Δ)
    end

    # check imag(z) << 1
    @assert maximum(imag(Ufilt)) < eps
    @assert minimum(imag(Ufilt)) > -eps
    @assert maximum(imag(Vfilt)) < eps
    @assert minimum(imag(Vfilt)) > -eps
    @assert maximum(imag(Wfilt)) < eps
    @assert minimum(imag(Wfilt)) > -eps

    return real(Ufilt), real(Vfilt), real(Wfilt)
end

function padvel(U::Array{Float32}, δ::Int64, n::Int64)
    np = 2*δ+n
    # Upad = Array{Float32}(undef, (2*δ+n, 2*δ+n, 2*δ+n))
    Upad = 12345*ones((np,np,np))
    for i in 1:np
        im = ind(i,n,δ)
        for j in 1:np
            jm = ind(j,n,δ)
            for k in 1:np
                km = ind(k,n,δ)
                Upad[i,j,k] = U[im,jm,km]
            end
        end
    end
    return Upad
end

function ind(i,n,δ)
    if i <= δ
        im = n - i
    elseif i > n
        im = i - n
    else 
        im = i
    end
    return im
end 

function totKE(U, V, W, n)::Float32
    KE = 0.
    for i in 1:n
        for j in 1:n
            for k in 1:n
                KE += 0.5*(U[i,j,k]^2 + V[i,j,k]^2 + W[i,j,k]^2)
            end
        end
    end
    return KE
end

function fluid2part(X,xm,ym,zm,U,V,W)::Vector{Float32}
    """
    Tri-linear velocity interpolation from Eulerian grid to Lagrangian particles
    """
    xp, yp, zp = X
    # Get the grid index
    dx = xm[2]-xm[1]
    dy = ym[2]-ym[1]
    dz = zm[2]-zm[1]
    nx = length(xm)
    ny = length(ym)
    nz = length(zm)
    ip=Int64(ceil(xp/dx))
    jp=Int64(ceil(yp/dy))
    kp=Int64(ceil(zp/dz))
    
    ip == 0 ? ip = 1 : nothing
    jp == 0 ? jp = 1 : nothing
    kp == 0 ? kp = 1 : nothing

    ip > nx ? ip = nx : nothing
    jp > ny ? jp = ny : nothing
    kp > nz ? kp = nz : nothing

    # Get interpolation cells w/r to centroids
    sx = -1
    sy = -1
    sz = -1
    if xp>=xm[ip]
        sx = 0
    end
    if yp>=ym[jp]
        sy = 0
    end
    if zp>=zm[kp]
        sz = 0
    end
    
    # Get the first interpolation point
    i1=ip+sx; i2=i1+1
    j1=jp+sy; j2=j1+1
    k1=kp+sz; k2=k1+1
    
    # Correct for periodicity
    if i1<1
        i1=nx
    end
    if i2>nx
        i2=1
    end
    if j1<1
        j1=ny
    end
    if j2>ny
        j2=1
    end
    if k1<1
        k1=nz
    end
    if k2>nz
        k2=1
    end
    
    # Compute the linear interpolation coefficients
    wx1 = (xm[i2]-xp )/(xm[i2]-xm[i1])
    wx2 = (xp -xm[i1])/(xm[i2]-xm[i1])
    wy1 = (ym[j2]-yp )/(ym[j2]-ym[j1])
    wy2 = (yp -ym[j1])/(ym[j2]-ym[j1])
    wz1 = (zm[k2]-zp )/(zm[k2]-zm[k1])
    wz2 = (zp -zm[k1])/(zm[k2]-zm[k1])

    wwd = zeros(2,2,2)
    
    # Combine the interpolation coefficients to form a tri-linear interpolation
    wwd[1,1,1]=wx1*wy1*wz1
    wwd[2,1,1]=wx2*wy1*wz1
    wwd[1,2,1]=wx1*wy2*wz1
    wwd[2,2,1]=wx2*wy2*wz1
    wwd[1,1,2]=wx1*wy1*wz2
    wwd[2,1,2]=wx2*wy1*wz2
    wwd[1,2,2]=wx1*wy2*wz2
    wwd[2,2,2]=wx2*wy2*wz2
    
    # Interpolate fluid velocity to particle location
    ug= wwd[1,1,1]*U[i1,j1,k1]+wwd[2,1,1]*U[i2,j1,k1]+wwd[1,2,1]*U[i1,j2,k1]+wwd[2,2,1]*U[i2,j2,k1]+wwd[1,1,2]*U[i1,j1,k2]+wwd[2,1,2]*U[i2,j1,k2]+wwd[1,2,2]*U[i1,j2,k2]+ wwd[2,2,2]*U[i2,j2,k2]
    vg= wwd[1,1,1]*V[i1,j1,k1]+wwd[2,1,1]*V[i2,j1,k1]+wwd[1,2,1]*V[i1,j2,k1]+wwd[2,2,1]*V[i2,j2,k1]+wwd[1,1,2]*V[i1,j1,k2]+wwd[2,1,2]*V[i2,j1,k2]+wwd[1,2,2]*V[i1,j2,k2]+ wwd[2,2,2]*V[i2,j2,k2]
    wg= wwd[1,1,1]*W[i1,j1,k1]+wwd[2,1,1]*W[i2,j1,k1]+wwd[1,2,1]*W[i1,j2,k1]+wwd[2,2,1]*W[i2,j2,k1]+wwd[1,1,2]*W[i1,j1,k2]+wwd[2,1,2]*W[i2,j1,k2]+wwd[1,2,2]*W[i1,j2,k2]+ wwd[2,2,2]*W[i2,j2,k2]
    return Vector{Float32}([ug, vg, wg])
end
