function createUpperLowerMapping(stgIndex::Int64, surf::TwoDSurfThick)

cam = zeros(length(surf.cam_slope[:]))
thick = zeros(length(surf.cam_slope[:]))
x = zeros(length(surf.cam_slope[:]))
qu = zeros(length(surf.cam_slope[:]))
ql = zeros(length(surf.cam_slope[:]))


cam[:] = surf.cam_slope[:]
thick[:] = surf.thick_slope[:]
x[:] = surf.x[:]
qu[:] = surf.ueU[:]
ql[:] = surf.ueL[:]


camber_slope_upper = [reverse(cam[2:stgIndex]);cam[2:end]]
thick_slope_upper = [reverse(thick[2:stgIndex]);thick[2:end]]

x_upper =  [-reverse(x[2:stgIndex]);x[1:end]]
x_lower =  x[stgIndex:end]

qustag =  [reverse(ql[2:stgIndex]);qu[2:end]]
qlstag  =  ql[stgIndex:end]

x_upper = x_upper .+ x[stgIndex]
x_lower = x_lower .- x[stgIndex]

camber_slope_lower = cam[stgIndex+1:end]
thick_slope_lower =  thick[stgIndex+1:end]

dsdx_upper = zeros(length(camber_slope_upper))
dsdx_lower = zeros(length(camber_slope_lower))


su = zeros(length(camber_slope_upper))
sl = zeros(length(camber_slope_lower))


for i = 1:length(camber_slope_upper)
        dsdx_upper[i] = sqrt(1 + (camber_slope_upper[i] + thick_slope_upper[i])^2)
end


for i = 1:length(camber_slope_lower)
        dsdx_lower[i] = sqrt(1 + (camber_slope_lower[i] + thick_slope_lower[i])^2)
end


su[1] = 0.
sl[1] = 0.

 for i = 1:length(camber_slope_upper)
       su[i] = simpleTrapz(dsdx_upper[1:i], x_upper[1:i])
 end

 for i = 1:length(camber_slope_lower)
       sl[i] = simpleTrapz(dsdx_lower[1:i], x_lower[1:i])
 end

return x_upper, x_lower, su, sl, qustag, qlstag

end
