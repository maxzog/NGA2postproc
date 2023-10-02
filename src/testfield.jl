
mutable struct testgrid<:grid
   n :: Int32
   L :: Float32
   Δ :: Float32
   U :: Array{Float32}
   V :: Array{Float32}
   W :: Array{Float32}
end

function init_testgrid(n::Int64, L::Float32)::testgrid
   delta = L / n
   xv = collect(LinRange(0.0f0, L - delta, n))
   
   funcU(x,y,z) = x^2 + y^2 + z^2
   funcV(x,y,z) = (x+3)^3 + sin(y) + sqrt(z)
   funcW(x,y,z) = x^4 + y^4 + z^4
   U = Array{Float32}(undef, (n, n, n))
   V = Array{Float32}(undef, (n, n, n))
   W = Array{Float32}(undef, (n, n, n))
   for i in 1:n, j in 1:n, k in 1:n
      U[i,j,k] = funcU(xv[i], xv[j], xv[k])
      V[i,j,k] = funcV(xv[i], xv[j], xv[k])
      W[i,j,k] = funcW(xv[i], xv[j], xv[k])
   end
   return testgrid(n, L, delta, U, V, W)
end

function eval_test(xi::Vector{Float32})
   funcU(x,y,z) = x^2 + y^2 + z^2
   funcV(x,y,z) = (x+3)^3 + sin(y) + sqrt(z)
   funcW(x,y,z) = x^4 + y^4 + z^4
   
   resU = funcU(xi[1], xi[2], xi[3])
   resV = funcV(xi[1], xi[2], xi[3])
   resW = funcW(xi[1], xi[2], xi[3])

   return [resU, resV, resW]
end
