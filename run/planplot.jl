include("../src/build.jl")

#using .Threads

function main(ARGS::Vector{String})
    dir = ARGS[1]
    ps = get_parts(dir*"/particles/", 90)
    plane_plot(ps, 90)
 #   @threads for i in 1:305
 #  	 ps  = get_parts(dir*"/particles/", i);
 #  	 #hit = get_grid(dir, i, 256, 6.2832)
 #  	 @printf("\nParticle file %i opened successfully", i)
 #  	 plane_plot(ps, i)
 #   end
end

function plane_plot(ps, i)
    # 2D slice
    ids     = slice_lagrangian(ps, "z", Float32(3.14), w=Float32(.1))
    x, y, z = getXs(ps[ids])
    #ωω = enstrophy(hit.U, hit.V, hit.W, hit.n, hit.n, hit.n, hit.L)
    #xm = collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2)
    #ym = collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2)
    #zm = collect(hit.Δ/2:hit.Δ:hit.L-hit.Δ/2)
    #println("Lastindex xm: ", lastindex(xm)) 
    #v = [scalar2part(p.pos, xm, ym, zm, ωω) for p in ps[ids]]
    writedlm("./Xs/x_"*string(i), x)
    writedlm("./Xs/y_"*string(i), y)
    writedlm("./Xs/z_"*string(i), z)
    #writedlm("./Xs/v_"*string(i), v)
    # plt = scatter(
    #         x, y, 
    #         markersize=.9, markercolor="black", α=1, 
    #         aspectratio=:equal, legend=:false, dpi=300,
    #         xlims=(0,6.2832), ylims=(0,6.2832), axis=([], false))
    # savefig("./figs/slice.png")
    return
end

function getXs(ps::Vector{part})
    np=lastindex(ps)
    x = zeros(Float32, np); y = zeros(Float32, np); z = zeros(Float32, np)
    for i in 1:lastindex(ps)
        x[i], y[i], z[i] = ps[i].pos
    end
    return x, y, z
end

function slice_lagrangian(ps::Vector{part}, plane::String, loc::Float32; w=Float32(0.1))
    ids = zeros(Int, lastindex(ps))
    for i in 1:lastindex(ps)
        p = ps[i]
        
        if plane=="x"
            if p.pos[1] < loc+w/2 && p.pos[1] > loc-w/2
                ids[i] = p.id
            end
        end
        if plane=="y"
            if p.pos[2] < loc+w/2 && p.pos[2] > loc-w/2
                ids[i] = p.id
            end
        end
        if plane=="z"
            if p.pos[3] < loc+w/2 && p.pos[3] > loc-w/2
                ids[i] = p.id
            end
        end
        
    end
    
    c = sum([id==0 ? 0 : 1 for id in ids])
    ids_slice = zeros(Int, c)
    
    ind = 1
    for i in 1:lastindex(ps)
        if ids[i] != 0
            ids_slice[ind] = ids[i]
            ind += 1
        end
    end
    return ids_slice
end

main(ARGS)
