#=
savefirst = "AOA_upper = "
plot(surf.x, surf.ueU)
savefig(string(savefirst,anglev))

close()
savefirst = "AOA_lower = "
plot(surf.x, surf.ueL)
savefig(string(savefirst,anglev))
close()

spuel = Spline1D(surf.x, surf.ueL)
print("lower rootsi :"roots(spuel)[1])
spueu = Spline1D(surf.x, surf.ueU)
print("upper roots : "roots(spueu)[1])
=#


B_mod = readdlm("modified_bspline.csv")
B_mod_1 = readdlm("modified_bspline_7thorder.csv")
B_ori = readdlm("original.csv")

ft = figure("thickness")
plot(B_ori[:,1], B_ori[:,2],"k", label="Original")
plot(B_mod[:,1], B_mod[:,2],"r--", label="Modified 6th order")
plot(B_mod_1[:,1], B_mod_1[:,2],"b--", label="Modified 7th order")
axis("equal")
legend()

fv = figure("edge velocity")
plot(B_ori[:,1], B_ori[:,3],"k")
plot(B_mod[:,1], B_mod[:,3],"r--")
plot(B_mod_1[:,1], B_mod_1[:,3],"r--")
a = fv.gca()
a.set_xlim([0.0,1.0])


fd = figure("delta")
plot(B_ori[:,1], B_ori[:,4],"k")
plot(B_mod[:,1], B_mod[:,4],"r--")
plot(B_mod_1[:,1], B_mod_1[:,4],"r--")

