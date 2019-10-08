push!(LOAD_PATH,"../UnsteadyFlowSolvers.jl/src/")
using UnsteadyFlowSolvers

cleanWrite()

#alphadef = EldRampReturnDef(45. *pi/180, 0.2, 11, 1.0)
alphadef = ConstDef(0. *pi/180)

hdef = ConstDef(0.0)

udef = ConstDef(1.)

full_kinem = KinemDef(alphadef, hdef, udef)

pvt = 0.0

geometry = "NACA0099"

lespcrit = [10.25;]

surf = TwoDSurfThick(geometry, pvt, full_kinem, ndiv=140, naterm=138)

curfield = TwoDFlowField()

dtstar = 0.015

t_tot = 7.0

nsteps = 1

println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = t_tot/5.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

#mat, surf, curfield = lautat_kutta(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)
mat, surf, curfield = lautat(surf, curfield, nsteps, dtstar,startflag, writeflag, writeInterval, delvort)

pop!(curfield.tev)

qu_base, ql_base = calc_edgeVel(surf, [0.; 0.])

i_sep = 47
xtev = surf.bnd_x_u[i_sep-3]+0.5*dt
ztev = surf.bnd_z_u[i_sep+3]
# ztevl = surf.bnd_z_l[surf.ndiv]
vcore = 0.02
tevstr = 0.05
tev1 = TwoDVort(xtev, ztev, tevstr, vcore, 0., 0.)
# tev2 = TwoDVort(xtev, ztevl, -tevstr, vcore, 0., 0.)
push!(curfield.tev, tev1)
# push!(curfield.tev, tev2)

update_indbound(surf, curfield)

qu_pair, ql_pair = calc_edgeVel(surf, [0.; 0.])

figure(1)
plot(surf.x, qu_base)
plot(surf.x, qu_pair)

figure(2)
plot(surf.x[2:end], diff(qu_base)./diff(surf.x))
plot(surf.x[2:end], diff(qu_pair)./diff(surf.x))
