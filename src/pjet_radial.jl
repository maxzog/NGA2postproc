
function get_radial_average(fn::String, step::Int64, nx::Int64, ny::Int64, nt::Int64)
   pjet = get_pjet(fn, step)
   # pjet.U .= 1.0f0
   # pjet.V .= 1.0f0
   # pjet.W .= 1.0f0
   xv = LinRange(0.0f0, 0.254f0, nx)
   rv = LinRange(-2.45f0 * 0.0254f0, 2.45f0 * 0.0254f0, ny)
   tmp = Array{Float32}(undef, (nx, ny))
   avg = radialgrid(nx, ny, maximum(xv), rv[end]-rv[1], xv, rv, 
                    Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)), 
                    Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)), 
                    Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)), Array{Float32}(undef, (nx, ny)))
   meanU = 0.0f0; meanV = 0.0f0; meanW = 0.0f0
   meanwX = 0.0f0; meanwY = 0.0f0; meanwZ = 0.0f0
   for i in 1:nx, j in 1:ny
      inds = get_radial_indices(pjet, xv[i], rv[j], nt)
      for k in 1:nt
         ii, jj, kk = inds[:,k]
         meanU += pjet.U[ii,jj,kk]
         meanV += pjet.V[ii,jj,kk]
         meanW += pjet.W[ii,jj,kk]
         meanwX += pjet.wX[ii,jj,kk]
         meanwY += pjet.wY[ii,jj,kk]
         meanwZ += pjet.wZ[ii,jj,kk]
      end
      avg.U[i,j] = meanU / nt
      avg.V[i,j] = meanV / nt
      avg.W[i,j] = meanW / nt
      avg.wX[i,j] = meanwX / nt
      avg.wY[i,j] = meanwY / nt
      avg.wZ[i,j] = meanwZ / nt
      meanU = 0.0f0; meanV = 0.0f0; meanW = 0.0f0
      meanwX = 0.0f0; meanwY = 0.0f0; meanwZ = 0.0f0
   end
   # meanU, meanV, meanW will now store RMS
   meanU = 0.0f0; meanV = 0.0f0; meanW = 0.0f0
   for i in 1:nx, j in 1:ny
      inds = get_radial_indices(pjet, xv[i], rv[j], nt)
      for k in 1:nt
         ii, jj, kk = inds[:,k]
         meanU += (pjet.U[ii,jj,kk]-avg.U[i,j])^2
         meanV += (pjet.V[ii,jj,kk]-avg.V[i,j])^2
         meanW += (pjet.W[ii,jj,kk]-avg.W[i,j])^2
      end
      avg.Urms[i,j] = sqrt(meanU / nt)
      avg.Vrms[i,j] = sqrt(meanV / nt)
      avg.Wrms[i,j] = sqrt(meanW / nt)
      meanU = 0.0f0; meanV = 0.0f0; meanW = 0.0f0
      meanwX = 0.0f0; meanwY = 0.0f0; meanwZ = 0.0f0
   end
   return avg
end

function get_lims(pjet::grid)
   ks = 1; ke = pjet.nz
   js = 1; je = pjet.ny
   is = 1; ie = pjet.nx
   
   isReplaced = false
   for i in 1:pjet.nx
      if pjet.xvm[i] > 0.0254f0*10.0f0 && !isReplaced
         ie = i
         isReplaced = true
      end
   end
   
   isReplaced = false
   isReplaced2 = false
   for i in 1:pjet.ny
      if pjet.yvm[i] > -0.0254f0*2.5f0 && !isReplaced
         js = i
         isReplaced = true
      end
      if pjet.yvm[i] > 0.0254f0*2.5f0 && !isReplaced2
         je = i
         isReplaced2 = true
      end
   end

   isReplaced = false
   isReplaced2 = false
   for i in 1:pjet.nz
      if pjet.zvm[i] > -0.0254f0*2.5f0 && !isReplaced
         ks = i
         isReplaced = true
      end
      if pjet.zvm[i] > 0.0254f0*2.5f0 && !isReplaced2
         ke = i
         isReplaced2 = true
      end
   end

   return [[is,ie], [js,je], [ks,ke]]
end

function get_radial_average(pjet::grid; nx_avg = 128, ny_avg = 64)
   avg = radialgrid(nx_avg, ny_avg, 10f0*0.0254f0, 2.5f0*0.0254f0, pjet.xvm, pjet.yvm, 
                    zeros(Float32, (nx_avg, ny_avg)), zeros(Float32, (nx_avg, ny_avg)), zeros(Float32, (nx_avg, ny_avg)), 
                    zeros(Float32, (nx_avg, ny_avg)), zeros(Float32, (nx_avg, ny_avg)), zeros(Float32, (nx_avg, ny_avg)), 
                    zeros(Float32, (nx_avg, ny_avg)))
   dr = maximum(pjet.yvm) / ny_avg
   dx = avg.Lx / avg.nx
   tmp = zeros(Int64, (nx_avg, ny_avg))
   lims = get_lims(pjet)
   is,ie = lims[1]
   js,je = lims[2]
   ks,ke = lims[3]
   for k in ks:ke
      for j in js:je
         for i in is:ie 
            r = [0.0f0, pjet.yvm[j], pjet.zvm[k]]
            rn = norm(r)
            rhat = r / rn
            rind = floor(Int64, rn / dr) + 1
            if !isapprox(norm(rhat), 1.0f0)
               rind = 1
            end
            rtheta = cross(rhat, [0.0f0,1.0f0,0.0f0])
            rtheta = cross(rtheta, rhat) / norm(cross(rtheta, rhat))
            omega = [pjet.wX[i,j,k], pjet.wY[i,j,k], pjet.wZ[i,j,k]]
            if !isapprox(norm(rtheta), 1.0f0)
               rtheta = omega / norm(omega)
            end

            i_avg = floor(Int64, pjet.xvm[i] / dx) + 1

            if rind < ny_avg && rind > 0 && i_avg < nx_avg
               avg.U[i_avg,rind] += pjet.U[i,j,k]
               avg.V[i_avg,rind] += pjet.V[i,j,k]
               avg.W[i_avg,rind] += pjet.W[i,j,k]
               avg.Urms[i_avg,rind] += pjet.U[i,j,k]^2
               avg.Vrms[i_avg,rind] += pjet.V[i,j,k]^2
               avg.Wrms[i_avg,rind] += pjet.W[i,j,k]^2
               avg.wT[i_avg,rind] += dot(omega, rtheta)
               tmp[i_avg,rind] += 1
            end
         end
         
      end
   end
   for j in 1:ny_avg
      for i in 1:nx_avg
         if tmp[i,j] != 0
            avg.U[i,j]  /= tmp[i,j]
            avg.V[i,j]  /= tmp[i,j]
            avg.W[i,j]  /= tmp[i,j]
            avg.Urms[i,j]  /= tmp[i,j]
            avg.Vrms[i,j]  /= tmp[i,j]
            avg.Wrms[i,j]  /= tmp[i,j]
            avg.wT[i,j] /= tmp[i,j]
         end
      end 
   end
   avg.Urms = avg.Urms - avg.U.^2
   avg.Vrms = avg.Vrms - avg.V.^2
   avg.Wrms = avg.Wrms - avg.W.^2
   return avg
end

function gen_points(x::Float32, r::Float32, ns::Int64)
   theta = LinRange(0.0, Float32(pi), ns)

   X = zeros(Float32, (ns,3))
   for i in eachindex(theta)
      X[i,1] = x
      X[i,2] = Float32(r*cos(theta[i]))
      X[i,3] = Float32(r*sin(theta[i]))
   end
   return X
end

function get_radial_indices(pjet::grid, x::Float32, r::Float32, ns::Int64)
   inds = Array{Int32}(undef, (3, ns))
   X = gen_points(x, r, ns)
   for s in 1:ns
      i,j,k = get_ijk_stretched(pjet, X[s, :])
      inds[1, s] = i
      inds[2, s] = j
      inds[3, s] = k
   end
   return inds
end

function get_ijk_stretched(pjet::grid, X::Vector{Float32}; maxiter = 10)
   i = index_inRange(X[1]+abs(pjet.xv[1]), pjet.dxm, pjet.xv .+ abs(pjet.xv[1]))
   j = index_inRange(X[2]+abs(pjet.yv[1]), pjet.dym, pjet.yv .+ abs(pjet.yv[1]))
   k = index_inRange(X[3]+abs(pjet.zv[1]), pjet.dzm, pjet.zv .+ abs(pjet.zv[1])) 
   return i, j, k
end

function index_inRange(x::Float32, dx::Float32, xv::Vector{Float32}; maxiter=25)::Int32
   i = floor(Int32, x / dx)
   if i == 0
      return 1
   end
   if i == 1
      i += 1
   end
   if i == lastindex(xv)
      i = lastindex(xv)-1
   end
   if i > lastindex(xv)
      return 99999999
   end
   if i < 0
      return -99999999
   end
   isIn = false; citer = 0
   while citer < maxiter && !isIn
      if i+1 > lastindex(xv)
         return i-1
      end
      if i-1 < 1
         return i+1
      end
      if x > xv[i-1] && x < xv[i+1]
         isIn = true
      elseif x < xv[i]
         i -= 1
      elseif x > xv[i]
         i +=1
      end
      citer += 1
   end
   return i
end

