push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers
using ForwardDiff
using PyPlot

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


phi_u = zeros(surf.ndiv)
phi_l = zeros(surf.ndiv)
cpu =   zeros(surf.ndiv)
cpl =  zeros(surf.ndiv)
t = 0.01
dt = 0.005
Re = 10000
cam_slope_orig = zeros(surf.ndiv)
isSep = false


pop!(curfield.tev)

#qu_base, ql_base = calc_edgeVel(surf, [0.; 0.])
qu_base, ql_base, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

i_sep = 47
xtev = surf.bnd_x_u[i_sep-3]+0.5*dt
ztev = surf.bnd_z_u[i_sep+3]
# ztevl = surf.bnd_z_l[surf.ndiv]
vcore = 0.02
tevstr = 0.002
tev1 = TwoDVort(xtev, ztev, tevstr, vcore, 0., 0.)
# tev2 = TwoDVort(xtev, ztevl, -tevstr, vcore, 0., 0.)
push!(curfield.tev, tev1)
# push!(curfield.tev, tev2)

update_indbound(surf, curfield)
qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

								       
#qu_pair, ql_pair = calc_edgeVel(surf, [0.; 0.])

smoothScaledEnd!(surf.x, qu,10)

thick_slope_orig = zeros(surf.ndiv)
delu, dell, Eu, El = initDelE(surf.ndiv-1)
delu_prev =  zeros(length(delu))
Eu_prev = zeros(length(delu))
thick_slope_orig[:] = surf.thick_slope[:] 
cam_slope_orig[:] = surf.cam_slope[:]
thick_orig = zeros(surf.ndiv) 
thick_orig[:] = surf.thick[:] 

s = zeros(surf.ndiv)

sc = zeros(surf.ndiv-1)

dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    s[1] = 0.
    for i = 2:surf.ndiv
        s[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
    end

 for i = 1:surf.ndiv-1
        sc[i] = (s[i] + s[i+1])/2
 end

quc = (qu[1:end-1] + qu[2:end])/2.0

qux = diff1(sc, quc)

qut = zeros(surf.ndiv-1)


wu = [delu delu.*(Eu .+ 1.0)] 
#error("stop")
wusoln, i_sepu, isSep = FVMIBLgridvar(wu, quc, qut, qux, diff(s), t-dt, t, isSep)

delu_prev[:] = delu[:]
delu[:] = wusoln[:,1]
Eu[:] = wusoln[:,2]./wusoln[:,1] .- 1.


smoothScaledEnd!(sc, delu, 10)


newthick = zeros(surf.ndiv)
con = zeros(surf.ndiv)
 
thickUpdate = zeros(surf.ndiv)

newcamb = zeros(surf.ndiv)
cambCon = zeros(surf.ndiv)
camUpdate = zeros(surf.ndiv)


 for i = 1:surf.ndiv-1

 con[i] =  (quc[i]*delu[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i] + cam_slope_orig[i]).^2))
             
 newthick[i] = thick_orig[i] + con[i]
  
  end

con[surf.ndiv] =  (quc[surf.ndiv-1]*delu[surf.ndiv-1])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[surf.ndiv] + cam_slope_orig[surf.ndiv]).^2))

 newthick[surf.ndiv] = thick_orig[surf.ndiv] + con[surf.ndiv]

bstart = [-0.1260; -0.3516; 0.2843; -0.1015; 0.0; 0.0]
coef = find_nacaCoef(surf, newthick, bstart)		   

th = parse(Int, surf.coord_file[7:8])/100.
b1 = 0.2969
b = [b1; coef]

@. nacath(x) = 5*th*(b[1]*sqrt(x) + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4 + b[6]*x^5 + b[7]*x^6)
 
 for i = 1:surf.ndiv
       surf.thick[i] = nacath(surf.x[i])
       surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
 end

qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

figure("Edge velocity")
plot(surf.x, qu_base, label= "Base")
plot(surf.x, qu, label = "Pair")
legend()

figure("Slope of velocity")
plot(surf.x[2:end], diff(qu_base)./diff(surf.x), label = "Base")
plot(surf.x[2:end], diff(qu)./diff(surf.x), label = "Pair")
legend()
