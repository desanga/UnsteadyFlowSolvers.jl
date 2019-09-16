
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
