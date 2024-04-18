
include("/Users/maxzog/codes/postproc/src/build.jl")

using .Threads

mutable struct Simulation
   # Time
   t     :: Float32
   dt    :: Float32
   tf    :: Float32
   step  :: Int
   # Domain
   L     :: Float32
   # Particles
   Npart :: Int
   Npair :: Int
end

global WIDTH = 0.04f0

function main()
   # Time
   t  = 0.0f0
   dt = Float32(1.0e-3)
   tf = 1.0f0
   # Domain
   L  = 3.14f0
   # Particles
   Npart = 100
   Npair = Int((Npart - 1) * Npart / 2)
   # Physics
   tau = 5.0f0
   sig = 0.5f0
   d   = 0.0f0
   taup= 0.1f0
  
   @threads for i in 1:10
      sim, ps = intialize(t, dt, tf, L, Npart, Npair, tau, sig, taup, d)
      
      @printf("\n-= Simulation parameters =-")
      @printf("\nNumber of particles - sim :: %i", sim.Npart)
      @printf("\nNumber of particles - ps  :: %i", lastindex(ps))
      @printf("\nSize of domain :: %.3f", sim.L)
      @printf("\ndt :: %.3f, tf :: %.3f", sim.t, sim.tf)
      @printf("\nTime :: %.3f", sim.t)
      @printf("\nStep :: %i", sim.step)
      
      sim, ps = pair_solve!(sim, ps)
      write_ps(ps, "./data/ou_parts_test_"*string(i)*".dat")
   end

end


function intialize(t, dt, tf, L, Npart, Npair, tau, sig, taup, d)
   sim = Simulation(t, dt, tf, 0, L, Npart, Npair)
   ps  = Vector{ou_part{Float32}}(undef, sim.Npart)
   for i in 1:sim.Npart
      pos = sim.L * rand(Float32, 3)
      ps[i] = ou_part(i, pos, [0.0f0,0.0f0,0.0f0],[0.0f0,0.0f0,0.0f0],[0.0f0,0.0f0,0.0f0],tau,sig,taup,d)
   end 
   update_dW!(ps)
   return sim, ps
end

function periodic!(p::ou_part, field::Simulation)
    for i in 1:3
        while p.pos[i] > field.L || p.pos[i] < 0.0
          if p.pos[i] > field.L
              p.pos[i] = p.pos[i] - field.L
          elseif p.pos[i] < 0.0
              p.pos[i] = p.pos[i] + field.L
          end
        end
    end
    return p
end

function S_f(r::Float32)
   t1 = WIDTH^3*pi^1.5*erf(40/WIDTH)^3
   t2 = 0.5*WIDTH^3*pi^1.5*erf(40/WIDTH)^2*(erf((40-r)/WIDTH)+erf((40+r)/WIDTH))
   t3 = 0.5*WIDTH^3*exp(-r^2/(2*WIDTH)^2)*pi^1.5*erf(40/WIDTH)^2*(erf((80-r)/2/WIDTH)+erf((80+r)/2/WIDTH))
   area = 2*sqrt(2)*WIDTH^3*pi^1.5*erf(20*sqrt(2)/WIDTH)^3
   sum = sqrt(2) * (t1 + t2 - 2*t3) / area
   return sum
end

function dS_f(r::Float32)
   t2 = 0.5*pi^1.5*WIDTH^3*erf(40/WIDTH)^2*(2*exp(-(r+40)^2/WIDTH^2)/sqrt(pi)/WIDTH - 2*exp(-(40-r)^2/WIDTH^2)/sqrt(pi)/WIDTH)
   t3 = 0.5*pi^1.5*WIDTH^3*erf(40/WIDTH)^2*exp(-r^2/4/WIDTH^2)*(exp(-(r+80)^2/4/WIDTH^2)/sqrt(pi)/WIDTH - exp(-(80-r)^2/4/WIDTH^2)/sqrt(pi)/WIDTH) - 
        0.25*pi^1.5*WIDTH*r*erf(40/WIDTH)^2*exp(-r^2/4/WIDTH^2)*(erf((80-r)/2/WIDTH)+erf((r+80)/WIDTH/2)) 
   area = 2*sqrt(2)*WIDTH^3*pi^1.5*erf(20*sqrt(2)/WIDTH)^3
   return sqrt(2)/area * (t2 - 2*t3)
end

function rho(p::ou_part{Float32}, q::ou_part{Float32})
   r  = q.pos - p.pos
   rn = norm(r)
   rh = r / rn
   return exp(-0.5 * rn^2 / WIDTH^2)
end

function B(p::ou_part{Float32}, q::ou_part{Float32})
   correlation_tensor = [1.0f0 0.0f0 0.0f0 rho(p,q) 0.0f0 0.0f0;
                         0.0f0 1.0f0 0.0f0 0.0f0 rho(p,q) 0.0f0;
                         0.0f0 0.0f0 1.0f0 0.0f0 0.0f0 rho(p,q);
                         rho(p,q) 0.0f0 0.0f0 1.0f0 0.0f0 0.0f0;
                         0.0f0 rho(p,q) 0.0f0 0.0f0 1.0f0 0.0f0;
                         0.0f0 0.0f0 rho(p,q) 0.0f0 0.0f0 1.0f0]
   return p.sig * sqrt(2*p.tau) * correlation_tensor * [p.dW; q.dW] / sqrt(1.0f0 + rho(p,q)^2)
end


function pair_solve!(sim::Simulation, ps::Vector{ou_part{Float32}})
   #=
   # This function solves many two-particle SCRW systems
   =#
   while sim.t < sim.tf
      for i in 1:2:lastindex(ps)
         p = ps[i]; q = ps[i+1] 
         rn = get_minr(p.pos, q.pos, sim.L)
         if rn < norm(q.pos-p.pos)
            rhat = -(q.pos - p.pos) / rn
         else
            rn = norm(q.pos - p.pos)
            rhat = (q.pos - p.pos) / rn
         end
         rn = norm(q.pos - p.pos)
         S_f_eval  = 2*p.sig^2 *  S_f(rn)
         dS_f_eval = 2*p.sig^2 * dS_f(rn)
        
         # Update fluid velocity seen
         p.fld += -p.fld*p.tau*sim.dt + B(p,q)[1:3]*sqrt(sim.dt) - (2*S_f_eval/rn + dS_f_eval)*(p.pos-q.pos)/rn*sim.dt 
         q.fld += -q.fld*q.tau*sim.dt + B(p,q)[4:6]*sqrt(sim.dt) + (2*S_f_eval/rn + dS_f_eval)*(q.pos-p.pos)/rn*sim.dt

         # Update particle velocity
         # -- check if particle is inertial or a tracer
         if p.d > 0.0f0
            p.vel += (p.fld - p.vel) / p.taup * sim.dt
            q.vel += (q.fld - q.vel) / q.taup * sim.dt
         else
            p.vel = p.fld
            q.vel = q.fld
         end

         # Update particle position
         p.pos += p.vel * sim.dt 
         q.pos += q.vel * sim.dt 
         periodic!(p, sim)
         periodic!(q, sim)
         ps[i] = p; ps[i+1] = q
      end

      update_dW!(ps)
      
      sim.t += sim.dt
      sim.step += 1

      print_state(sim, ps)
   end
   return sim, ps
end

function print_state(sim::Simulation, ps::Vector{ou_part{Float32}})
   @printf("\nTime :: %.3f ----- Step :: %5i", sim.t, sim.step)
end

function update_dW!(ps::Vector{ou_part{Float32}})
   for i in 1:lastindex(ps)
      p = ps[i]
      p.dW = randn(Float32, 3)
      ps[i] = p
   end
end

# function solve!(sim::Simulation, ps::Vector{ou_part{Float32}})
#    #=
#    # This function solves many two-particle SCRW systems
#    =#
#    while sim.t < sim.tf
#       for i in 1:lastindex(ps)
#          p = ps[i]
#          for j in i+1:lastindex(ps)
#             q = ps[j] 
#             rn = get_minr(p.pos, q.pos, sim.L)
#             if rn < norm(q.pos-p.pos)
#                rhat = -(q.pos - p.pos) / rn
#             else
#                rhat = (q.pos - p.pos) / rn
#             end
#             S_f_eval  = 2*p.sig^2 *  S_f(rn)
#             dS_f_eval = 2*p.sig^2 * dS_f(rn)
#             
#             # Update fluid velocity seen
#             # del_r \cdot <delta u_i delta u_j> = (2 * S_f(r) / r + d S_f(r) / dr) rhat
#             # Process:
#             # 1. Evaluate structure function and derivative at r
#             # 2. Dot with separation vector
#             p.fld += B(p,q)[1:3]*sqrt(sim.dt) + (2*S_f_eval/rn + dS_f(rn))*rhat*sim.dt 
#          end
#          p.fld += -p.fld*p.tau*sim.dt
#          # Update particle velocity
#          # -- check if particle is inertial or a tracer
#          if p.d > 0.0f0
#             p.vel += (p.fld - p.vel) / p.taup * sim.dt
#          else
#             p.vel = p.fld
#          end
#
#          # Update particle position
#          p.pos += p.vel * sim.dt 
#          periodic!(p, sim)
#          ps[i] = p
#       end
#       
#       sim.t += sim.dt
#       sim.step += 1
#
#       print_state(sim, ps)
#    end
#    return sim, ps
# end

@time main()
