function reconstructGrid(stgIndex::Int64, surf::TwoDSurfThick, su::Array{Float64}, qu::Array{Float64}, ql::Array{Float64}, qu_prev::Array{Float64}, ql_prev::Array{Float64})

x = zeros(length(surf.cam_slope[:]))


x[:] = surf.x[:]

x_u =  [-reverse(x[2:stgIndex]);x[1:end]]
x_l =  x[stgIndex:end]

qustag =  [reverse(ql[2:stgIndex]);qu[1:end]]
qlstag  =  ql[stgIndex:end]


qustag_prev =  [reverse(ql_prev[2:stgIndex]);qu_prev[1:end]]
qlstag_prev  =  ql_prev[stgIndex:end]

x_u = x_u .+ x[stgIndex]
x_l = x_l .- x[stgIndex]

s_u =  [-reverse(su[2:stgIndex]);su[1:end]]
s_l =  su[stgIndex:end]

s_u = s_u .+ su[stgIndex]
s_l = s_l .- su[stgIndex]

return x_u, x_l, s_u, s_l, qustag, qlstag, qustag_prev, qlstag_prev

end

function reverseReconstructGrid(stgIndex::Int64, delu::Array{Float64}, dell::Array{Float64}, qustag::Array{Float64}, qlstag::Array{Float64})

	deluu = zeros(surf.ndiv-1)
	delll = zeros(surf.ndiv-1)
	quc = zeros(surf.ndiv-1)
	qlc = zeros(surf.ndiv-1)


	deluu[:] = delu[stgIndex:end];
	delll[:] = [delu[1:stgIndex];dell[1:end]]

	qu[:] = qustag[stgIndex:end];
	ql[:] = [qustag[1:stgIndex];qlstag[1:end]]




	return deluu, delll, quc, qlc
end

