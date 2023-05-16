
function w_r(ps::Vector{part}, nb::Int64; L=Float32(6.2832))
    wr = zeros(Float32, (2, nb))
    cs = zeros(Int, (2, nb))
    rmax = norm([L/2, L/2, L/2])
    dr = rmax/nb
    for i in 1:lastindex(ps)
        for j in i+1:lastindex(ps)
            r̂ = unit_r(ps[i].pos, ps[j].pos)
            w = dot(ps[j].vel-ps[i].vel, r̂)
            r = get_minr(ps[j].pos, ps[i].pos, L)
            #r = norm(ps[j].pos - ps[i].pos)
	    rind = floor(Int, r/dr) + 1
            if w > 0
                wr[1,rind] += w
                cs[1,rind] += 1
            elseif w < 0
                wr[2,rind] += w
                cs[2,rind] += 1
            end
        end
    end
#    for i in 1:nb
#        if cs[1,i] != 0
#            wr[1,i] /= cs[1,i]
#        end
#        if cs[2,i] != 0
#            wr[2,i] /= cs[2,i]
#        end
#    end
    return wr, cs
end

function s_r(ps::Vector{part}, nb::Int64; L=Float32(6.2832))
    wr = zeros(Float32, nb)
    cs = zeros(Int, nb)
    rmax = norm([L/2, L/2, L/2])
    dr = rmax/nb
    for i in 1:lastindex(ps)
        for j in i+1:lastindex(ps)
            r̂ = unit_r(ps[i].pos, ps[j].pos)
            w = dot(ps[j].fld-ps[i].fld, r̂)
            r = get_minr(ps[j].pos, ps[i].pos, L)
	        rind = floor(Int, r/dr) + 1
            wr[rind] += w*w
            cs[rind] += 1
        end
    end
#    for i in 1:nb
#        if cs[1,i] != 0
#            wr[1,i] /= cs[1,i]
#        end
#        if cs[2,i] != 0
#            wr[2,i] /= cs[2,i]
#        end
#    end
    return wr, cs
end

function s_r_op(ps::Vector{part}, nb::Int64; L=Float32(6.2832))
    wr = zeros(Float32, nb)
    cs = zeros(Int, nb)
    rmax = norm([L/2, L/2, L/2])
    dr = rmax/nb
    for i in 1:1
        for j in i+1:lastindex(ps)
            r̂ = unit_r(ps[i].pos, ps[j].pos)
            w = dot(ps[j].fld-ps[i].fld, r̂)
            r = get_minr(ps[j].pos, ps[i].pos, L)
	        rind = floor(Int, r/dr) + 1
            wr[rind] += w*w
            cs[rind] += 1
        end
    end
    return wr, cs
end

function unit_r(x1::Vector{Float32}, x2::Vector{Float32})::Vector{Float32}
    return (x2-x1)/norm(x2-x1)
end

