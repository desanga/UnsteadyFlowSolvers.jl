push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers
using DelimitedFiles

cleanWrite()

alphadef = ConstDef(5.0*pi/180)

hdef = ConstDef(0.)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.25

#geometry = "Cylinder"
geometry = "NACA0012"

lespcrit = [10.25;]

surf = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=136)

curfield = TwoDFlowField()

dtstar = 0.005
t_tot = 2.5

nsteps = Int(round(t_tot/dtstar))+1
println("nsteps ", nsteps)

startflag = 0

writeflag = 0

writeInterval = 0.1

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

Re = 10000

mat, surf, curfield =  IBL_shape_attached(Re, surf, curfield, nsteps, dtstar, startflag, writeflag, writeInterval, delvort)



#mat, surf, curfield = lautat(surf,curfield, nsteps ,dtstar, startflag, writeflag, writeInterval,delvort, maxwrite = 100, nround=6)

