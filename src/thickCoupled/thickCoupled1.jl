function createUpperLowerMapping(stgIndex::Int64, surf::TwoDSurfThick)

cam = zeros(length(surf.cam_slope[:]))
thick = zeros(length(surf.cam_slope[:]))
x = zeros(length(surf.cam_slope[:]))
qu = zeros(length(surf.cam_slope[:]))
ql = zeros(length(surf.cam_slope[:]))
#camber_slope = zeros(stgIndex)
#thick_slope = zeros(stgIndex)


cam[:] = surf.cam[:]
thick[:] = surf.thick[:]
x[:] = surf.x[:]
qu[:] = surf.ueU[:]
ql[:] = surf.ueL[:]

x_u =  [-reverse(x[2:stgIndex]);x[1:end]]
x_l =  x[stgIndex:end]

qustag =  [reverse(ql[2:stgIndex]);qu[1:end]]
qlstag  =  ql[stgIndex:end]

x_u = x_u .+ x[stgIndex]
x_l = x_l .- x[stgIndex]

camber_stag_lower = cam[stgIndex:end]
thick_stag_lower =  thick[stgIndex:end]
camber_stag_upper = [reverse(cam[2:stgIndex]);cam[1:end]]
thick_stag_upper = [reverse(thick[2:stgIndex]);thick[1:end]]


camber_slope_upper = diff(camber_stag_upper)./diff(x_u)
thick_slope_upper = diff(thick_stag_upper)./diff(x_u)

camber_slope_lower = diff(camber_stag_lower)./diff(x_l)
thick_slope_lower = diff(thick_stag_lower)./diff(x_l)


dsdx_upper = zeros(length(camber_slope_upper))
dsdx_lower = zeros(length(camber_slope_lower))


su = zeros(length(x_u))
sl = zeros(length(x_l))


for i = 1:length(camber_slope_upper)
        dsdx_upper[i] = sqrt(1 + (camber_slope_upper[i] + thick_slope_upper[i])^2)
end


for i = 1:length(camber_slope_lower)
        dsdx_lower[i] = sqrt(1 + (camber_slope_lower[i] + thick_slope_lower[i])^2)
end


su[1] = 0.
sl[1] = 0.

 for i = 2:length(camber_slope_upper)
       su[i] = simpleTrapz(dsdx_upper[1:i], x_u[1:i])
 end

 for i = 2:length(camber_slope_lower)
       sl[i] = simpleTrapz(dsdx_lower[1:i], x_l[1:i])
 end

return x_u, x_l, su, sl, qustag, qlstag

end

function createUpperLowerRemapping(stgIndex::Int64, delu::Array{Float64}, dell::Array{Float64})





return delu, dell, qu, ql
end

