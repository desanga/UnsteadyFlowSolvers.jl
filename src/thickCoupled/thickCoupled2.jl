function reconstructGrid(stgIndex::Int64, surf::TwoDSurfThick, su::Array{Float64}, qu::Array{Float64}, ql::Array{Float64}, qu_prev::Array{Float64}, ql_prev::Array{Float64}, delu::Array{Float64}, dell::Array{Float64}, Eu::Array{Float64}, El::Array{Float64}, delu_iter::Array{Float64}, dell_iter::Array{Float64}, delu_prev::Array{Float64}, dell_prev::Array{Float64}, Eu_iter::Array{Float64}, El_iter::Array{Float64})

x = zeros(length(surf.cam_slope[:]))
x[:] = surf.x[:]

x_u =  [-reverse(x[2:stgIndex]);x[1:end]]
x_l =  x[stgIndex:end]
x_u = x_u .+ x[stgIndex]
x_l = x_l .- x[stgIndex]

s_u =  [-reverse(su[2:stgIndex]);su[1:end]]
s_l =  su[stgIndex:end]
s_u = s_u .+ su[stgIndex]
s_l = s_l .- su[stgIndex]

qustag =  [reverse(ql[2:stgIndex]);qu[1:end]]
qlstag  =  ql[stgIndex:end]

qustag_prev =  [reverse(ql_prev[2:stgIndex]);qu_prev[1:end]]
qlstag_prev  =  ql_prev[stgIndex:end]


delustag = [reverse(dell[2:stgIndex]);delu[1:end]]
Eustag =  [reverse(El[2:stgIndex]);Eu[1:end]]

dellstag  =  dell[stgIndex:end]
Elstag  =  El[stgIndex:end]

delu_iterstag = [reverse(dell_iter[2:stgIndex]);delu_iter[1:end]]
dell_iterstag = dell_iter[stgIndex:end]

delu_prevstag = [reverse(dell_prev[2:stgIndex]);delu_prev[1:end]]
dell_prevstag = dell_prev[stgIndex:end]

Eu_iterstag = [reverse(El_iter[2:stgIndex]);Eu_iter[1:end]]
El_iterstag =  El_iter[stgIndex:end]


return x_u, x_l, s_u, s_l, qustag, qlstag, qustag_prev, qlstag_prev, delustag, dellstag, Eustag, Elstag, dell_iterstag, delu_iterstag, delu_prevstag, dell_prevstag, Eu_iterstag, El_iterstag

end

function reverseReconstructGrid(stgIndex::Int64, surf::TwoDSurfThick, qustag::Array{Float64}, qlstag::Array{Float64}, qustagc::Array{Float64}, qlstagc::Array{Float64}, qustagc_prev::Array{Float64}, qlstagc_prev::Array{Float64}, delustag::Array{Float64}, dellstag::Array{Float64}, Eustag::Array{Float64}, Elstag::Array{Float64}, dellstag_iter::Array{Float64}, delustag_iter::Array{Float64}, delustag_prev::Array{Float64}, dellstag_prev::Array{Float64}, Eustag_iter::Array{Float64}, Elstag_iter::Array{Float64})


	quc = zeros(surf.ndiv-1)
	qlc = zeros(surf.ndiv-1)
	delu = zeros(surf.ndiv-1)
	dell = zeros(surf.ndiv-1)
	Eu = zeros(surf.ndiv-1)
	El = zeros(surf.ndiv-1)
	quc_prev = zeros(surf.ndiv-1)
	qlc_prev = zeros(surf.ndiv-1)
	quc_prev = zeros(surf.ndiv-1)
	qlc_prev = zeros(surf.ndiv-1)
	qu  =  zeros(surf.ndiv)
	ql = zeros(surf.ndiv)

	quc_prev =  zeros(surf.ndiv-1)
	qlc_prev =  zeros(surf.ndiv-1)


	delu_iter =  zeros(surf.ndiv-1)
	dell_iter =  zeros(surf.ndiv-1)

	delu_prev =  zeros(surf.ndiv-1)
	dell_prev =  zeros(surf.ndiv-1)

	Eu_iter =  zeros(surf.ndiv-1)
	El_iter =  zeros(surf.ndiv-1)


	delu[:] = delustag[stgIndex:end]
	dell[:] = [reverse(delustag[1:stgIndex]);dellstag[2:end]]

	Eu[:] = Eustag[stgIndex:end]
	El[:] = [reverse(Eustag[1:stgIndex]);Elstag[2:end]]


	quc[:] = qustagc[stgIndex:end]
	qlc[:] = [reverse(qustagc[1:stgIndex]);qlstagc[2:end]]

	 
	qu[:] = qustag[stgIndex:end]
	ql[:] = [reverse(qustag[1:stgIndex]);qlstagc[1:end]]
	 
	quc_prev[:] = qustagc_prev[stgIndex:end]
	qlc_prev[:] = [reverse(qlstagc_prev[1:stgIndex]);qlstagc_prev[2:end]]
	
	delu_iter[:] = delustag_iter[stgIndex:end]
	dell_iter[:] = [reverse(delustag_iter[1:stgIndex]);dellstag_iter[2:end]] 
	
	delu_prev[:] = delustag_prev[stgIndex:end] 
	dell_prev[:] = [reverse(delustag_prev[1:stgIndex]);dellstag_prev[2:end]] 
	
	Eu_iter[:] = Eustag_iter[stgIndex:end]
	El_iter[:] = [reverse(Eustag_iter[1:stgIndex]);Elstag_iter[2:end]]


	return quc, qlc, quc_prev, qlc_prev, delu, dell, Eu, El, delu_iter, dell_iter, delu_prev, dell_prev, Eu_iter, El_iter

	end


