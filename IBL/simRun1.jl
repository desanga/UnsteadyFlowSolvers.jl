push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers
using DelimitedFiles


Re = parse(Float64,ARGS[1])
geometry = ARGS[2]
tot = parse(Float64,ARGS[3])
output = ARGS[4]

cleanWrite()

alphadef = ConstDef(0. *pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

#geometry = "Cylinder"
#geometry = "NACA0012"

lespcrit = [10.25;]

surf = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=136)

curfield = TwoDFlowField()

dtstar = 0.005
t_tot = tot

nsteps = Int(round(t_tot/dtstar))+1
println("nsteps ", nsteps)

startflag = 0

writeflag = 0

writeInterval = 0.1

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

#Re = 10000

mat, surf, curfield = IBL_shape_attached(Re, surf, curfield, nsteps, dtstar, startflag, writeflag, writeInterval, delvort)

delta_scaled = (surf.delta/sqrt(Re))
ue = surf.Ue
xx = surf.x
thick = surf.thick
data = [xx ue delta_scaled thick]
open(output,"w") do io
	writedlm(io,data)
end
