include("build.jl")

function main()
    # dir = "/Users/maxzog/NGA2/examples/hit/ensight/HIT/particles/"
    dir = "/Users/maxzog/cases/hit_les/ensight/HIT/particles/"
    vt = 0.; c = 0
    for i in 150:200
        ps = get_parts(dir, i)
        vt += variance([p.fld for p in ps])
        c  += 1
    end
    @printf("Variance = %.3f", vt/c)
    return
end


main()