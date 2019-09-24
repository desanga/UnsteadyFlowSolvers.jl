camber_slope_upper = [reverse(surf.cam_slope[2:5]);surf.cam_slope[2:end]]
thick_slope_upper = [reverse(surf.thick_slope[2:5]);surf.thick_slope[2:end]]

x_upper =  [-reverse(surf.x[2:5]);surf.x[2:end]]
x_lower =  surf.x[6:end]

x_upper = x_upper .+ surf.x[5]
x_lower = x_lower .- surf.x[6]


println(length(camber_slope_upper))
println(length(camber_slope_lower))

camber_slope_lower = surf.cam_slope[6:end]
thick_slope_lower =  surf.thick_slope[6:end]

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






