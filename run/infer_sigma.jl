include("../src/build.jl")

function main()
	psn = get_parts("/home/maxzog/NGA2s/nga2restart/examples/hit/ensight/HIT/particles/", 5)
	psm = get_parts("/home/maxzog/NGA2s/nga2restart/examples/hit/ensight/HIT/particles/", 4)
	hit = get_grid("/home/maxzog/NGA2s/nga2restart/examples/hit/ensight/HIT", 5, 256, 6.2832)
	
	ids = sample(1:1_000_000, 350_000, replace=false, ordered=true)
	@time dr, sig, s = pic_cov(psn[ids], psm[ids], hit, 1.2, 256)
	writedlm("./sigma_ar.txt",sig)
	writedlm("./sigma_sf.txt",  s)
	writedlm("./sigma_dr.txt", dr)
end


main()
