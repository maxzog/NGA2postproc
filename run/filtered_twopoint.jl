
include("../src/build.jl")

function main()
    dir = "/home/maxzog/NGA2/examples/hit/ensight_tracer/HIT/"
    step=13
    hit = get_grid(dir, step, 256, 6.2832)
    deltas = [4, 8, 16, 32]

    ns=1_000; rmax=3.14; nbins=128
    for d in deltas
        les = LES(hit, "sharpspec", d)
        psd = get_parts_dns(dir*"particles/", 13)
        psf = get_parts_dns(dir*"particles/", 13)
        pss = get_parts_dns(dir*"particles/", 13)
     
        Usgs = les.U-les.Uf
        Vsgs = les.V-les.Vf
        Wsgs = les.W-les.Wf

        xm = collect(les.Δ/2:les.Δ:les.L-les.Δ/2)

        for i in 1:lastindex(psd)
            p = psf[i]; q = pss[i]
            uf=fluid2part(p.pos, xm, xm, xm, les.Uf, les.Vf, les.Wf)
            us=fluid2part(q.pos, xm, xm, xm, Usgs, Vsgs, Wsgs)
            p.fld=Float32.(uf); q.fld=Float32.(us)
            psf[i]=p; pss[i]=q
        end

        dr, dlf, dtf, sf = pic_struct_lt(sample(psf, ns, replace=false), les, rmax, bins)
        dr, dls, dts, ss = pic_struct_lt(sample(pss, ns, replace=false), les, rmax, bins)

        writedlm("./outs/structs/filtered_"*string(d)*"_t.dat", dtf)
        writedlm("./outs/structs/filtered_"*string(d)*"_l.dat", dlf)
     
        writedlm("./outs/structs/sgs_"*string(d)*"_t.dat", dts)
        writedlm("./outs/structs/sgs_"*string(d)*"_l.dat", dls)
     end
end

main()
