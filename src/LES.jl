
function spectral!(hit::lesgrid)
    p = plan_fft(hit.U)
    hit.Uk = p*hit.U; hit.Vk = p*hit.V; hit.Wk = p*hit.W; hit.Pk = p*hit.P
    return 
end

function ispectral!(hit::lesgrid)
    ip = plan_ifft(hit.Uf)
    hit.Uk=ip*hit.Uk; hit.Vk=ip*hit.Vk; hit.Wk=ip*hit.Wk; hit.Pk=ip*hit.Pk
    return 
end

function tophat_spec!(hit::lesgrid, kv::Vector{Int64})
    for i∈1:hit.n, j∈1:hit.n, k∈1:hit.n
        kk = norm([kv[i] kv[j] kv[k]])
        if kk > 1e-6
            G = sin(0.5*kk*hit.Δf)/(0.5*kk*hit.Δf)
            hit.Uk[i,j,k] *= G 
            hit.Vk[i,j,k] *= G 
            hit.Wk[i,j,k] *= G 
            hit.Pk[i,j,k] *= G
        end
    end
    return 
end

function sharpspec!(hit::lesgrid, kv::Vector{Int64})
    for i∈1:hit.n, j∈1:hit.n, k∈1:hit.n
        kk = norm([kv[i] kv[j] kv[k]])
        heavy = hit.κc - kk
        if heavy <= 0
            hit.Uk[i,j,k] = 0.0
            hit.Vk[i,j,k] = 0.0
            hit.Wk[i,j,k] = 0.0
            hit.Pk[i,j,k] = 0.0
        end
    end
    return 
end

function gaussian_spec!(hit::lesgrid, kv::Vector{Int64})
    for i∈1:hit.n, j∈1:hit.n, k∈1:hit.n
        kk = norm([kv[i] kv[j] kv[k]])
        G  = exp(-kk^2*hit.Δf^2/24)
        hit.Uk[i,j,k] *= G
        hit.Vk[i,j,k] *= G
        hit.Wk[i,j,k] *= G
        hit.Pk[i,j,k] *= G
    end
    return 
end

function spectralfilter!(hit::lesgrid)
    spectral!(hit)
    n = trunc(Int, hit.n/2)
    kv = fftshift(-n:n-1)
    if hit.ftype == "sharpspec"
       sharpspec!(hit, kv)
    elseif hit.ftype == "gaussian"
       gaussian_spec!(hit, kv) 
    elseif hit.ftype == "tophat"
       tophat_spec!(hit, kv) 
    end
    ispectral!(hit)
    hit.isFiltered=true
    return 
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

function LES(hit::dnsgrid, ftype::String, δ::Int64; eps = 1e-5)
   Δf = δ * hit.L / hit.n
   les = lesgrid(ftype, false, hit.n, hit.L, hit.Δ, Δf, π/Δf, hit.U, hit.V, hit.W, hit.P,
                 zeros(Float32, hit.n, hit.n, hit.n),zeros(Float32, hit.n, hit.n, hit.n),
                 zeros(Float32, hit.n, hit.n, hit.n),zeros(Float32, hit.n, hit.n, hit.n),
                 zeros(Complex{Float32}, hit.n, hit.n, hit.n),zeros(Complex{Float32}, hit.n, hit.n, hit.n),
                 zeros(Complex{Float32}, hit.n, hit.n, hit.n),zeros(Complex{Float32}, hit.n, hit.n, hit.n))
   spectralfilter!(les)
   # check imag(z) << 1
   println(maximum(imag(les.Uk)))
   @assert maximum(imag(les.Uk)) < eps
   @assert minimum(imag(les.Uk)) > -eps
   @assert maximum(imag(les.Vk)) < eps
   @assert minimum(imag(les.Vk)) > -eps
   @assert maximum(imag(les.Wk)) < eps
   @assert minimum(imag(les.Wk)) > -eps

   les.Uf=real(les.Uk); les.Vf=real(les.Vk); les.Wf=real(les.Wk); les.Pf=real(les.Pk);
   return les
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
    return KE/(n^3)
end

# THIS DOES NOT WORK
function data2part(X,hit::grid)::Vector{Float64}
    vel = Vector{Float64}(undef, 3)

    ip=ceil(Int64, X[1]/hit.Δ)
    jp=ceil(Int64, X[2]/hit.Δ) 
    kp=ceil(Int64, X[3]/hit.Δ)

    x=collect(0:hit.n)*hit.Δ
    y=collect(0:hit.n)*hit.Δ
    z=collect(0:hit.n)*hit.Δ
    
    xm=collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2) 
    ym=collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2) 
    zm=collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2) 
    
    # Interpolate U
    i = maximum(minimum(hit.n-1,ip), 1)
    while X[1]-x[i]<0 && i>1
        i-=1
    end
    while X[1]-x[i+1]>=0 && i+1<hit.n
        i+=1
    end
    
    j = maximum(minimum(hit.n-1,jp), 1)
    while X[2]-ym[j]<0 && j>1
        j-=1
    end
    while X[2]-ym[j+1]>=0 && j+1<hit.n
        j+=1
    end

    k = maximum(minimum(hit.n-1,kp), 1)
    while X[3]-zm[k]<0 && k>1
        k-=1
    end
    while X[3]-zm[k+1]>=0 && k+1<hit.n
        k+=1
    end
    wx1=(X[1]-x[ip] )/hit.Δ; wx2=1.0-wx1
    wy1=(X[2]-ym[ip])/hit.Δ; wy2=1.0-wy1
    wz1=(X[3]-zm[ip])/hit.Δ; wz2=1.0-wz1

    vel[1]=wz1*(wy1*(wx1*hit.U[ip+1,jp+1,kp+1] + wx2*hit.U[ip  ,jp+1,kp+1])  + 
                wy2*(wx1*hit.U[ip+1,jp  ,kp+1] + wx2*hit.U[ip  ,jp  ,kp+1])) + 
           wz2*(wy1*(wx1*hit.U[ip+1,jp+1,kp  ] + wx2*hit.U[ip  ,jp+1,kp  ])  + 
                wy2*(wx1*hit.U[ip+1,jp  ,kp  ] + wx2*hit.U[ip  ,jp  ,kp  ]))
    
    # Interpolate V
    wx1=(X[1]-xm[ip])/hit.Δ; wx2=1.0-wx1
    wy1=(X[2]-y[ip] )/hit.Δ; wy2=1.0-wy1
    wz1=(X[3]-zm[ip])/hit.Δ; wz2=1.0-wz1

    vel[2]=wz1*(wy1*(wx1*hit.V[ip+1,jp+1,kp+1] + wx2*hit.V[ip  ,jp+1,kp+1])  + 
                wy2*(wx1*hit.V[ip+1,jp  ,kp+1] + wx2*hit.V[ip  ,jp  ,kp+1])) + 
           wz2*(wy1*(wx1*hit.V[ip+1,jp+1,kp  ] + wx2*hit.V[ip  ,jp+1,kp  ])  + 
                wy2*(wx1*hit.V[ip+1,jp  ,kp  ] + wx2*hit.V[ip  ,jp  ,kp  ]))
    
    # Interpolate W
    wx1=(X[1]-xm[ip])/hit.Δ; wx2=1.0-wx1
    wy1=(X[2]-ym[ip])/hit.Δ; wy2=1.0-wy1
    wz1=(X[3]-z[ip] )/hit.Δ; wz2=1.0-wz1

    vel[3]=wz1*(wy1*(wx1*hit.W[ip+1,jp+1,kp+1] + wx2*hit.W[ip  ,jp+1,kp+1])  + 
                wy2*(wx1*hit.W[ip+1,jp  ,kp+1] + wx2*hit.W[ip  ,jp  ,kp+1])) + 
           wz2*(wy1*(wx1*hit.W[ip+1,jp+1,kp  ] + wx2*hit.W[ip  ,jp+1,kp  ])  + 
                wy2*(wx1*hit.W[ip+1,jp  ,kp  ] + wx2*hit.W[ip  ,jp  ,kp  ]))

    return vel
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

function scalar2part(X,xm,ym,zm,U)::Float32
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
    return ug
end
