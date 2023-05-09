#! /usr/bin/env/ julia

include("../src/build.jl")

using Plots

mutable struct simulation
    Δt        :: Float32
    t         :: Float32
    tf        :: Float32
    eps       :: Float32
    τLi       :: Float32
    τpi       :: Float32
    corr_func :: String
end

function main()
    dt=1e-3; tf=2.0; func="Hello"; 
    np=5_000; L=6.2832; 
    eps=27.0; tke=15.0; nu=8e-3; 
    St=1
    sim, hit, ps = init_sim(dt, tf, func, np, L, eps, tke, nu, St)
    step = 0

    make_gif = true

    @printf("\nStarting simulation ...")
    while sim.t < sim.tf
        # sim_step_scrw!(sim, ps, hit)
        sim_step_rel!(sim, ps, hit)
        sim.t += sim.Δt
        @printf("\nTime = %.3f", sim.t)
        if step%5==0
            save_parts(sim, ps)
        end
        step += 1
    end
    @printf("\nSimulation finished")
    if make_gif
        a  = Animation()
        fs = glob("./outs/langevin/*")
        for i in 1:lastindex(fs)
            X = readdlm(fs[i])
            frame(a, part_scatter(X, sim, hit))
        end
        gif(a, "./img/animation.gif")
    end
end

function part_scatter(X, sim::simulation, hit::grid)
    plt = scatter(X[:,1], X[:,2], X[:,3], 
                  xlim=(0,hit.L), ylim=(0,hit.L), zlim=(0,hit.L), legend=false,
                  markersize=0.5, markerstrokewidth=.5, fα=0.8, mα=0.8, dpi=300, aspect_ratio=:equal)
    return plt
end

function save_parts(sim::simulation, ps::Vector{part})
    writedlm("./outs/langevin/pos_"*string(round(sim.t, digits=3)), [p.pos for p in ps])
end

function sim_step_rel!(sim::simulation, ps::Vector{part}, field::grid)
    p = ps[1]
    p.pos = Float32.([field.L/2, field.L/2, field.L/2])
    p.vel = Float32.([0., 0., 0.])
    p.fld = Float32.([0., 0., 0.])
    p.uf  = Float32.([0., 0., 0.])
    ps[1] = p
    for i in 2:lastindex(ps)
        q = ps[i] 
        dW = randn(Float32, 6)
       
        rm  = get_minr(p.pos,ps[i].pos,field.L)
        rho = corr(p, ps[i], sim, field)
        
        B   = [1-rho 0 0 rho-1 0 0 ;
               0 1-rho 0 0 rho-1 0 ;
               0 0 1-rho 0 0 rho-1]
        
        q.fld = q.fld - q.fld*sim.τLi*sim.Δt + sqrt(1.5*sim.eps)*B*dW*sqrt(sim.Δt)/sqrt(1+rho^2)
        q.vel = q.vel + (q.fld - q.vel)*sim.τpi*sim.Δt
        q.pos = q.pos + q.vel*sim.Δt

        periodic!(q, field)

        ps[i] = q
    end
end

function sim_step_scrw!(sim::simulation, ps::Vector{part}, field::grid)
    for i in 1:lastindex(ps)
        # Spatial correlation part - particle-in-cell?
        # Could just sum over all particles
        
        p = ps[i]
        BdWi=zeros(Float32, 3); den=0.0; dSdr=zeros(Float32, 3)
       
        for j in 1:lastindex(ps)
            BdWi += sqrt(1.5*sim.eps*sim.Δt)*ps[j].uf*corr(p, ps[j], sim, field)
            r  = ps[j].pos-p.pos
            rm = get_minr(p.pos,ps[j].pos,field.L)
            if i != j
                dSdr += 2/3*sim.eps^(2/3)*rm^(-1/3)*r/norm(r)
            end
            den  += corr(p, ps[j], sim, field)^2
        end

        # dS/dr = 2/3*eps^(2/3)*r^(-1/3)*r/norm(r)
        
        p.fld = p.fld - p.fld*sim.τLi*sim.Δt + dSdr*sim.Δt/lastindex(ps) + BdWi/sqrt(den)
        p.vel = p.vel + (p.fld - p.vel)*sim.τpi*sim.Δt
        p.pos = p.pos + p.vel*sim.Δt

        periodic!(p, field)

        ps[i] = p
    end
    update_dW!(ps)
end

function corr(p::part, q::part, sim::simulation, field::grid)
    r = get_minr(p.pos, q.pos, field.L)
    r_c = 0.17; sig = 2*.1^2
    rhoij = (exp(-r^2/sig) - exp(-r_c^2/sig)) / (1-exp(r_c^2/sig))
    return rhoij
end

function periodic!(p::part, field::grid)
    for i in 1:3
        if p.pos[i] > field.L
            p.pos[i] = p.pos[i] - field.L
        elseif p.pos[i] > field.L
            p.pos[i] = p.pos[i] + field.L
        end
    end
    return p
end

function rand_pos!(ps::Vector{part}, field::grid)
    for i in 1:lastindex(ps)
        ps[i].pos = rand(Float32, 3)*field.L
    end
    return ps
end

function update_dW!(ps::Vector{part})
    for i in 1:lastindex(ps)
        ps[i].uf = randn(Float32, 3)
    end
    return ps
end

function init_sim(dt::Float64, tf::Float64, func::String, np::Int64, L::Float64, 
                  eps::Float64, tke::Float64, nu::Float64, St::Float64)
    τLi = (0.5+0.75*2.1)*eps/tke
    dp  = sqrt(18*nu*sqrt(nu/eps)*St/1000)
    τp  = 1000*dp^2/18/nu; τpi = 1/τp
    sim = simulation(Float32(dt), Float32(0), Float32(tf), Float32(eps), Float32(τLi), Float32(τpi), func)
    ps  = init_parts(np)
    update_dW!(ps)
    field = dnsgrid(64, Float32(L), Float32(L/64), zeros(Float32, (3,3,3)), 
                    zeros(Float32, (3,3,3)), zeros(Float32, (3,3,3)), zeros(Float32, (3,3,3)))
    rand_pos!(ps, field)
    return sim, field, ps
end

main()
