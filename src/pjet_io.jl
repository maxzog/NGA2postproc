
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
    xvm :: Vector{Float32}
    yvm :: Vector{Float32}
    zvm :: Vector{Float32}
    dxm :: Float32
    dym :: Float32
    dzm :: Float32
    U :: Array{Float32}
    V :: Array{Float32}
    W :: Array{Float32}
    wX :: Array{Float32}
    wY :: Array{Float32}
    wZ :: Array{Float32}
end

mutable struct pipegrid<:grid
    nx :: Int64
    ny :: Int64
    nz :: Int64
    Lx :: Float32
    Ly :: Float32
    Lz :: Float32
    D  :: Float32
    xv :: Vector{Float32}
    yv :: Vector{Float32}
    zv :: Vector{Float32}
    dxm :: Float32
    dym :: Float32
    dzm :: Float32
    U :: Array{Float32}
    V :: Array{Float32}
    W :: Array{Float32}
end

mutable struct radialgrid<:grid
   nx :: Int64
   ny :: Int64
   Lx :: Float32
   Ly :: Float32
   xv :: Vector{Float32}
   yv :: Vector{Float32}
   U :: Array{Float32}
   V :: Array{Float32}
   W :: Array{Float32}
   wT :: Array{Float32}
   Urms :: Array{Float32}
   Vrms :: Array{Float32}
   Wrms :: Array{Float32}
end

function get_pjet(fn::String, step::Int64)
   """
   Reads NGA2 ensight files and returns a grid type containing the stored data fields
   """

   step = string(step)

   nx, ny, nz, xv, yv, zv = read_mesh(fn)
   Lx = xv[end]-xv[1]
   Ly = yv[end]-yv[1]
   Lz = zv[end]-zv[1]

   xvm = (xv[2:end] .+ xv[1:end-1]) / 2
   yvm = (yv[2:end] .+ yv[1:end-1]) / 2
   zvm = (zv[2:end] .+ zv[1:end-1]) / 2
   

   dxm = mean(xv[2:end] - xv[1:end-1])
   dym = mean(yv[2:end] - yv[1:end-1])
   dzm = mean(zv[2:end] - zv[1:end-1])

   fnv = fn*"/velocity/velocity." * "0"^(6-length(step)) * step
   fnw = fn*"/vorticity/vorticity." * "0"^(6-length(step)) * step

   U, V, W = read_vel(fnv, nx, ny, nz)
   vortX, vortY, vortZ = read_vel(fnw, nx, ny, nz)
   return pjetgrid(nx, ny, nz, Lx, Ly, Lz, xv, yv, zv, xvm, yvm, zvm, dxm, dym, dzm,
                   U, V, W, vortX, vortY, vortZ)
end

function get_pipe(fn::String, step::Int64) 
   """
   Reads NGA2 ensight files and returns a grid type containing the stored data fields
   """
   
   step = string(step)
   nx, ny, nz, xv, yv, zv = read_mesh(fn)
   Lx = xv[end]-xv[1]
   Ly = yv[end]-yv[1]
   Lz = zv[end]-zv[1]
   
   dxm = xv[2] - xv[1]
   dym = yv[2] - yv[1]
   dzm = zv[2] - zv[1]

   fnv = fn*"/velocity/velocity." * "0"^(6-length(step)) * step
   U, V, W = read_vel(fnv, nx, ny, nz)

   return pipegrid(nx, ny, nz, Lx, Ly, Lz, 0.0f0, xv, yv, zv, dxm, dym, dzm,
                   U, V, W)
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

function write_avg(avg::radialgrid, step::Int64; dir = "./avg_out/")
   suf = string(step)*".txt"
   writedlm(dir*"U_"*suf, avg.U)
   writedlm(dir*"V_"*suf, avg.V)
   writedlm(dir*"W_"*suf, avg.W)
   writedlm(dir*"Urms_"*suf, avg.Urms)
   writedlm(dir*"Vrms_"*suf, avg.Vrms)
   writedlm(dir*"Wrms_"*suf, avg.Wrms)
   writedlm(dir*"wT_"*suf, avg.wT)
   writedlm(dir*"Lx_"*suf, avg.Lx)
   writedlm(dir*"Ly_"*suf, avg.Ly)
   return 
end
