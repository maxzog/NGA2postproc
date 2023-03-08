
include("../src/build.jl")

using .Threads

function par_wr()
	ps = get_parts("/home/maxzog/NGA2s/nga2restart/examples/hit/ensight/HIT/particles/", 80)
	nb = 128
	wp = zeros(Float32, (nthreads(), 2, nb)) 
	cp = zeros(Int64,   (nthreads(), 2, nb))
	@threads for i in 1:25
		psp = sample(ps, 25_000, replace=false, ordered=true)
		w, c = w_r(psp, nb)
		wp[threadid(), :, :] += w
		cp[threadid(), :, :] += c
	end
	w = reshape(sum(wp, dims=1), 2, nb)
	c = reshape(sum(cp, dims=1), 2, nb)
	writedlm("./relvel_dns_test.txt", w)
	writedlm("./relvel_dns_test_c.txt", c)
end

@time par_wr()

