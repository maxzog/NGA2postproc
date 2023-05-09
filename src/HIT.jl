function ∇(arr, n, Δxi, Δyi, Δzi, inds)
    """
    Second-order accurate central difference
    """
    i, j, k = inds[1]
    ip1, jp1, kp1 = inds[2]
    im1, jm1, km1 = inds[2]
    ddx = (arr[ip1, j, k] - 2*arr[i,j,k] + arr[im1, j, k])*Δxi/2
    ddy = (arr[i, jp1, k] - 2*arr[i,j,k] + arr[i, jm1, k])*Δyi/2
    ddz = (arr[i, j, kp1] - 2*arr[i,j,k] + arr[i, j, km1])*Δzi/2
    return ddx, ddy, ddz
end

function enstrophy(U::Array{Float32, 3}, V::Array{Float32, 3}, W::Array{Float32, 3},
    nx::Int64, ny::Int64, nz::Int64, L::AbstractFloat)
    """
    Returns enstrophy field
    Assumes triply periodic cubic domain
    """
    ωω = zeros(nx, ny, nz)
    dxi = nx/L; dyi = ny/L; dzi = nz/L
    for i ∈ 1:nx, j ∈ 1:ny, k ∈ 1:nz
        ip1 = i + 1; im1 = i - 1
        jp1 = j + 1; jm1 = j - 1
        kp1 = k + 1; km1 = k - 1
        i == nx ? ip1 = 1 : nothing
        j == ny ? jp1 = 1 : nothing
        k == nz ? kp1 = 1 : nothing
        i == 1 ? im1 = nx : nothing
        j == 1 ? jm1 = ny : nothing
        k == 1 ? km1 = nz : nothing

        inds = [(i,j,k), (ip1,jp1,kp1), (im1,jm1,km1)]

        dudx, dudy, dudz = ∇(U, nx, dxi, dyi, dzi, inds)
        dvdx, dvdy, dvdz = ∇(V, ny, dxi, dyi, dzi, inds)
        dwdx, dwdy, dwdz = ∇(W, nz, dxi, dyi, dzi, inds)
        ω = [dwdy-dvdz, dudz-dwdx, dvdx-dudy]
        ωω[i,j,k] = dot(ω, ω)
    end
    return ωω
end

function dissipation(field::grid, nu::Float64, veltype::String)
    """
    Returns psuedo-dissipation rate
    """
    if veltype == "SGS"
        Usgs = field.U - field.Uf
        Vsgs = field.V - field.Vf
        Wsgs = field.W - field.Wf
    end
    eps_p = 0
    dxi = 1/field.Δ; dyi = 1/field.Δ; dzi = 1/field.Δ;
    for i ∈ 1:field.n, j ∈ 1:field.n, k ∈ 1:field.n
        ip1 = i + 1; im1 = i - 1
        jp1 = j + 1; jm1 = j - 1
        kp1 = k + 1; km1 = k - 1
        i == field.n ? ip1 = 1 : nothing
        j == field.n ? jp1 = 1 : nothing
        k == field.n ? kp1 = 1 : nothing
        i == 1 ? im1 = field.n : nothing
        j == 1 ? jm1 = field.n : nothing
        k == 1 ? km1 = field.n : nothing

        inds = [(i,j,k), (ip1,jp1,kp1), (im1,jm1,km1)]
        if veltype == "DNS"
            dudx, dudy, dudz = ∇(field.U, field.n, dxi, dyi, dzi, inds)
            dvdx, dvdy, dvdz = ∇(field.V, field.n, dxi, dyi, dzi, inds)
            dwdx, dwdy, dwdz = ∇(field.W, field.n, dxi, dyi, dzi, inds)
        elseif veltype == "SGS"
            dudx, dudy, dudz = ∇(Usgs, field.n, dxi, dyi, dzi, inds)
            dvdx, dvdy, dvdz = ∇(Vsgs, field.n, dxi, dyi, dzi, inds)
            dwdx, dwdy, dwdz = ∇(Wsgs, field.n, dxi, dyi, dzi, inds)
        elseif veltype == "LES"
            dudx, dudy, dudz = ∇(field.Uf, field.n, dxi, dyi, dzi, inds)
            dvdx, dvdy, dvdz = ∇(field.Vf, field.n, dxi, dyi, dzi, inds)
            dwdx, dwdy, dwdz = ∇(field.Wf, field.n, dxi, dyi, dzi, inds)
        end
        eps_p += dudx^2 + dudy^2 + dudz^2 + 
                 dvdx^2 + dvdy^2 + dvdz^2 + 
                 dwdz^2 + dwdy^2 + dwdz^2
    end
    return nu*eps_p/field.n^3
end
