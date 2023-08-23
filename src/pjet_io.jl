
mutable struct pjetgrid<:grid
    nx :: Int64
    ny :: Int64
    nz :: Int64
    Lx :: Float32
    Ly :: Float32
    Lz :: Float32
    xv :: Vector{Float32}
    yv :: Vector{Float32}
    zv :: Vector{Float32}
    U :: Array{Float32}
    V :: Array{Float32}
    W :: Array{Float32}
    P :: Array{Float32}
end

function get_pjet(fn::String, step::Int64)
   """
   Reads NGA2 binary data files and returns a grid type containing the stored data fields
   """

   step = string(step)

   nx, ny, nz, xv, yv, zv = read_mesh(fn)
   Lx = xv[end]-xv[1]
   Ly = yv[end]-yv[1]
   Lz = zv[end]-zv[1]

   fnv = fn*"/velocity/velocity." * "0"^(6-length(step)) * step
   fnp = fn*"/pressure/pressure." * "0"^(6-length(step)) * step

   U, V, W = read_vel(fnv, nx, ny, nz)
   # P = read_P(fnp, nx, ny, nz)
   P = Array{Float32}(undef, (nx, ny, nz))
   return pjetgrid(nx, ny, nz, Lx, Ly, Lz, xv, yv, zv, 
                   U, V, W, P)
end
   

function read_mesh(fn::String)
   ## Read the mesh
   io = open(fn*"/geometry", "r")
   # skip ensight crap
   cbuff = Vector{Char}(undef, 80)
   for _ in 1:6
      for i in eachindex(cbuff)
         cbuff[i] = read(io, Char)
      end
   end
   # skip extents
   read!(io, Vector{Float32}(undef, 6))
   # skip more ensight crap
   for i in eachindex(cbuff)
      cbuff[i] = read(io, Char)
   end
   read!(io, Vector{Int32}(undef, 1))
   for _ in 1:2 
      for i in eachindex(cbuff)
         cbuff[i] = read(io, Char)
      end
   end

   # read mesh
   ncells = Vector{Int32}(undef, 3)
   read!(io, ncells)
   xv = Vector{Float32}(undef, ncells[1])
   yv = Vector{Float32}(undef, ncells[2])
   zv = Vector{Float32}(undef, ncells[3])
   read!(io, xv)
   read!(io, yv)
   read!(io, zv)
   close(io)
   nx = Int32(ncells[1]-1)
   ny = Int32(ncells[2]-1)
   nz = Int32(ncells[3]-1)
   return nx, ny, nz, xv, yv, zv 
end
