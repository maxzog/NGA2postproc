
function get_jhgrid(dir::String, n::Int64, L::Float64)
   tmpU = Array{Float32}(undef, (n,n,n))
   tmpV = Array{Float32}(undef, (n,n,n))
   tmpW = Array{Float32}(undef, (n,n,n))

   io = open(dir*"Um.dat", "r")
   read!(io, tmpU); close(io) 
   
   io = open(dir*"Vm.dat", "r")
   read!(io, tmpV); close(io) 
   
   io = open(dir*"Wm.dat", "r")
   read!(io, tmpW); close(io) 

   hit = dnsgrid(n, Float32(L), Float32(L/n), tmpU, tmpV, tmpW, Array{Float32}(undef, (n,n,n)))
   return hit
end
