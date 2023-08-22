
function get_minr(X1::Vector{Float32}, X2::Vector{Float32}, perL::Float32; eps = 1e-10)::Float32
    """
    Minimum particle separation distance assuming cubic, periodic domain
    """
    r = zeros(length(X1))
    for i ∈ 1:length(X1)
        dx = X2[i] - X1[i]
        
        # abs(dx) > abs(dx - perL) ? dx -= perL : nothing
        # abs(dx) > abs(dx + perL) ? dx += perL : nothing
        r[i] = dx - perL*round(Int, dx/perL)
    end
    return norm(r)
end

function spectral(hit::grid)
    p = plan_fft(hit.U)
    Uk=p*hit.U; Vk=p*hit.V; Wk=p*hit.W
    return Uk, Vk, Wk 
end

function spectral(hit::lesgrid)
    p = plan_fft(hit.Uf)
    Uk=p*hit.Uf; Vk=p*hit.Vf; Wk=p*hit.Wf
    return Uk, Vk, Wk 
end

function ispectral(Uk::Array{Float32}, Vk::Array{Float32}, Wk::Array{Float32})
    ip = plan_ifft(Uk)
    U=ip*Uk; V=ip*Vk; W=ip*Wk
    return U, V, W 
end

function get_spectrum(hit::grid)  
    n = trunc(Int, hit.n/2)
    kv = fftshift(-n:n-1)
    dk = 2*π/hit.L
    Uk, Vk, Wk = spectral(hit)
    spec = zeros(Float32, ceil(Int64, norm([maximum(abs.(kv))+.5 for _ in 1:3])))
    counts = zeros(Int64, length(spec))
    for i∈1:hit.n, j∈1:hit.n, k∈1:hit.n
        kk = norm([kv[i], kv[j], kv[k]])
        ik = floor(Int64, kk/dk+0.5)+1
        if ik > 1e-6
            spec[ik] += 0.5*(real(Uk[i,j,k]*conj(Uk[i,j,k])) + real(Vk[i,j,k]*conj(Vk[i,j,k])) + real(Wk[i,j,k]*conj(Wk[i,j,k])))
        end
    end
    return spec/hit.n^6
end

function get_spectrum(hit::lesgrid)  
    n = trunc(Int, hit.n/2)
    kv = fftshift(-n:n-1)
    dk = 2*π/hit.L
    Uk, Vk, Wk = spectral(hit)
    spec = zeros(Float32, ceil(Int64, norm([maximum(abs.(kv))+.5 for _ in 1:3])))
    counts = zeros(Int64, length(spec))
    for i∈1:hit.n, j∈1:hit.n, k∈1:hit.n
        kk = norm([kv[i], kv[j], kv[k]])
        ik = floor(Int64, kk/dk+0.5)+1
        if ik > 1e-6
            spec[ik] += 0.5*(real(Uk[i,j,k]*conj(Uk[i,j,k])) + real(Vk[i,j,k]*conj(Vk[i,j,k])) + real(Wk[i,j,k]*conj(Wk[i,j,k])))
        end
    end
    return spec/hit.n^6
end

function get_tauEM(hitn::grid, hitm::grid, dt::Float32)::Float32
   num = mean(hitn.U .* hitn.U + 
              hitn.V .* hitn.V + 
              hitn.W .* hitn.W)
   
   Ax = (hitn.U - hitm.U) ./ dt
   Ay = (hitn.V - hitm.V) ./ dt
   Az = (hitn.W - hitm.W) ./ dt
   
   den = mean(Ax .* Ax + 
              Ay .* Ay +
              Az .* Az)
   return sqrt(num/den)
end

function get_tauEM(hitnp1::grid, hitn::grid, hitnm1::grid, dt::Float32)::Float32
   num = mean(hitn.U .* hitn.U + 
              hitn.V .* hitn.V + 
              hitn.W .* hitn.W)
   
   Ax = (hitnp1.U - hitnm1.U) ./ (2 * dt)
   Ay = (hitnp1.V - hitnm1.V) ./ (2 * dt)
   Az = (hitnp1.W - hitnm1.W) ./ (2 * dt)
   
   den = mean(Ax .* Ax + 
              Ay .* Ay +
              Az .* Az)
   return sqrt(num/den)
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

function Lag4Interp(field::grid, xiv::Vector{Float32})::Vector{Float32}   
   n = floor(xiv[1] / field.Δ)
   p = floor(xiv[2] / field.Δ)
   q = floor(xiv[3] / field.Δ)
   xv = LinRange(0.0f0, field.L, field.n)
   rvec = zeros(Float32, 3)
   for i in 1:4, j in 1:4, k in 1:4
      ii = n-2+i
      jj = p-2+j
      kk = q-2+k
      lx = get_Lag4poly(xv, n, xiv[1], ii)
      ly = get_Lag4poly(xv, p, xiv[2], jj)
      lz = get_Lag4poly(xv, q, xiv[3], kk)
      rvec[1] += field.U[ii, jj, kk] * lx * ly * lz
      rvec[2] += field.V[ii, jj, kk] * lx * ly * lz
      rvec[3] += field.W[ii, jj, kk] * lx * ly * lz
   end
   return rvec
end

function get_Lag4poly(xv::Vector{Float32}, ind::Int64, xi::Float32, i::Int64)
   num = 1.0; den = 1.0
   for j in ind-1:ind+2
      if i != j
         num *= (xi - xv[j])
         den *= (xv[i] - xv[j])
      end
   end
   return num/den
end   
