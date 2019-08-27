push!(LOAD_PATH,"../src/")
using UnsteadyFlowSolvers

using LsqFit
using ForwardDiff
using Dierckx
using Serialization
using DelimitedFiles
using PyPlot
## Add any new functions here



alphadef = ConstDef(0. *pi/180)

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

t_tot = 5.

nsteps = Int(round(t_tot/dtstar))+1
println("nsteps ", nsteps)

startflag = 0

writeflag = 1

writeInterval = 0.5

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()


### Code below should go inside the solver function

function IBLnew(surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)
    
    # If a restart directory is provided, read in the simulation data
    if startflag == 0
        mat = zeros(0, 12)
        t = 0.
    elseif startflag == 1
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
        latestTime = maximum(dirresults)
        mat = readdlm("resultsSummary")
        t = mat[end,1]
    else
        throw("invalid start flag, should be 0 or 1")
    end
    mat = mat'

    dt = dtstar*surf.c/surf.uref

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        tTot = nsteps*dt
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dt)))
            end
        end
    end

    vcore = 1.3*dt*surf.c

    int_wax = zeros(surf.ndiv)
    int_c = zeros(surf.ndiv)
    int_t = zeros(surf.ndiv)

    RHStransp = zeros(surf.ndiv*2-2)

    wtu = zeros(surf.ndiv)
    wtl = zeros(surf.ndiv)

    t = 0.
    dt = dtstar

    Re = 10000

    #At this time coded for symmetric situations only
    qu = zeros(surf.ndiv)
    su = zeros(surf.ndiv)
    quc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)
    qucx = zeros(surf.ndiv-1)
    quct = zeros(surf.ndiv-1)
    suc = zeros(surf.ndiv-1)

    thick_orig = zeros(surf.ndiv)
    thick_orig[:] = surf.thick[:]
    thick_slope_orig = zeros(surf.ndiv)
    thick_slope_orig[:] = surf.thick_slope[:]

    #Initialise boundary layer
    delu, _, Eu, _ = initDelE(surf.ndiv-1)

    dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = 2:surf.ndiv
        su[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
   end

    for i = 1:surf.ndiv-1
        suc[i] = (su[i] + su[i+1])/2
    end

    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)

    #wk_src = TwoDSource[]

    x_w = collect(surf.c*1.01:surf.c*0.01:surf.c*3)
    nw = length(x_w)
    wfn = zeros(nw)

    fit_d = zeros(5)

    bl_end = surf.ndiv-1

    iterMax = 20
    
    for istep = 1:nsteps

        #Udpate current time
        t = t + dt

        #Update kinematic parameters
        update_kinem(surf, t)

        #Update flow field parameters if any
        update_externalvel(curfield, t)

        #Update bound vortex positions
        update_boundpos(surf, dt)

        #Update incduced velocities on airfoil
        update_indbound(surf, curfield)


        iter = 0
        res = 1
        del_iter = zeros(surf.ndiv-1)
        del_prev = zeros(surf.ndiv-1)
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(surf.ndiv-1)

        resTol = 1e-5

        quc_prev[:] = quc[:]       

        i_sep = 0
        
        while res > resTol
            
            iter += 1

            #Set up the matrix problem
            surf, xloc_tev, zloc_tev = update_thickLHS(surf, curfield, dt, vcore)
            
            #Construct RHS vector
            update_thickRHS(surf, curfield)
            
            
            soln = surf.LHS[1:surf.ndiv*2-3, 1:surf.naterm*2+1] \ surf.RHS[1:surf.ndiv*2-3]
            
            #Assign the solution
            for i = 1:surf.naterm
                surf.aterm[i] = soln[i]
                surf.bterm[i] = soln[i+surf.naterm]
            end
            tevstr = soln[2*surf.naterm+1]*surf.uref*surf.c
            
            temptev = TwoDVort(xloc_tev, zloc_tev, soln[2*surf.naterm+1], vcore, 0., 0.)

            #ind_new_u_u, ind_new_w_u = ind_vel_src(wk_src, surf.bnd_x_u, surf.bnd_z_u)
            #ind_new_u_l, ind_new_w_l = ind_vel_src(wk_src, surf.bnd_x_l, surf.bnd_z_l)
            
            #surf.uind_u[:] += ind_new_u_u[:]
            #surf.wind_u[:] += ind_new_w_u[:]
            #surf.uind_l[:] += ind_new_u_l[:]
            #surf.wind_l[:] += ind_new_w_l[:]
          
            a_wfn = log10(2*quc[end]*del_iter[end]/sqrt(Re))
            wfn[1] = 2*(quc[end]*del_iter[end] - quc[end-1]*del_iter[end-1])/(sqrt(Re)*(surf.x[end] - surf.x[end-1]))
            for i = 2:nw
                wfn[i] = (10^(a_wfn - 3.2*(x_w[i] - surf.c)) - 10^(a_wfn - 3.2*(x_w[i] - 1.)))/(x_w[i] - x_w[i-1])
            end
            
            #Source strengths in wake and induced velocity
            uind_src = zeros(surf.ndiv)
            for i = 1:surf.ndiv
                for iw = 1:nw-1
                    str = 0.5(wfn[iw] + wfn[iw+1])
                    xloc = 0.5*(x_w[iw] + x_w[iw+1])
                    uind_src[i] += 1/(2*pi)*str/(xloc - surf.x[i])*(x_w[iw+1] - x_w[iw])
                end
                # for iw = bl_end:surf.ndiv-2
                #     str = (quc[iw+1]*del_iter[iw+1] - quc[iw]*del_iter[iw])/(sqrt(Re)*(surf.x[iw+1] - surf.x[iw]))
                #     xloc = 0.5*(surf.x[iw] + surf.x[iw+1])
                #     uind_src[i] += 1/(2*pi)*str/(xloc - surf.x[i])*(x_w[iw+1] - x_w[iw])
                # end
            end
            
            surf.uind_u[:] += uind_src[:]
            surf.uind_l[:] += uind_src[:]
            
            ind_new_u_u, ind_new_w_u = ind_vel([temptev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]

            qu[:], _ = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

            surf.uind_u[:] -= uind_src[:]
            surf.uind_l[:] -= uind_src[:]
            
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]

            #smoothEnd!(qu,10)
            smoothScaledEnd!(surf.x,qu,10)
            
            #Solve the FV problem at cell centres
            for i = 1:surf.ndiv-1
                quc[i] = (qu[i] + qu[i+1])/2
            end
            
            qucx[2:end] = diff(quc)./diff(suc)
            qucx[1] = qucx[2]

            #smoothEnd!(qucx, 10)
            
            if istep == 1
                quct[:] .= 0.
            else
                quct[:] = (quc[:] - quc_prev[:])/dt
            end

            w = [delu delu.*(Eu .+ 1)]
          
            wsoln, i_sep = FVMIBLgridvar(w[1:bl_end,:], quc[1:bl_end], quct[1:bl_end], qucx[1:bl_end], diff(su)[1:bl_end], t-dt, t) 
            
            del_prev[:] = del_iter[:]
            
            del_iter[1:bl_end] = wsoln[:,1]
           
            E_iter[1:bl_end] = wsoln[:,2]./wsoln[:,1] .- 1.

            # if bl_end != surf.ndiv-1
            #     del_iter[bl_end+1:end] .= del_iter[bl_end] 
            #     E_iter[bl_end+1:end] .= E_iter[bl_end]
            # end
            
            if i_sep != 0
                println("separation detected")
                #bl_end = i_sep-1
                #smoothScaledEnd!(suc, del_iter, surf.ndiv-i_sep-2)
                #smoothScaledEnd!(suc, E_iter, surf.ndiv-i_sep-2)
            end
            
            #smoothEnd!(del_iter, 10)
            smoothScaledEnd!(suc, del_iter, 10)
            
            #Find suitable naca coefficients to fit the modified airfoil

            newthick = zeros(surf.ndiv)
            for i = 1:surf.ndiv-1
                newthick[i] = thick_orig[i] + (quc[i]*del_iter[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i]).^2))
            end
                        
            newthick[surf.ndiv] = thick_orig[surf.ndiv] + (quc[surf.ndiv-1]*del_iter[surf.ndiv-1])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[surf.ndiv]).^2))

            bstart = [-0.1260; -0.3516; 0.2843; -0.1015]
           
            coef = find_nacaCoef(surf, newthick, bstart)
            
            th = parse(Int, surf.coord_file[7:8])/100.
            b1 = 0.2969
            b = [b1; coef]

            @. nacath(x) = 5*th*(b[1]*sqrt(x) + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4)
            
            #Find new shape of airfoil
            for i = 1:surf.ndiv
                surf.thick[i] = nacath(surf.x[i])
                surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
            end
            
            #Check for convergence
            res =  sum(abs.(del_prev .- del_iter))
            println(iter, "   ", res)

            #if iter == iterMax
            if res <= resTol
                println("converged")
                
                # #Re-mesh the delta to prevent spurious oscillations developing
                # @. deltacurve(x, d) = d[1] + d[2]*x + d[3]*x^2 + d[4]*x^3 + d[5]*x^4
                # fit = curve_fit(deltacurve, suc[:], del_iter[:], fit_d)
                # fit_d = coef(fit)
                # @. deltacurve(x) = fit_d[1] + fit_d[2]*x + fit_d[3]*x^2 + fit_d[4]*x^3 + fit_d[5]*x^4
                # for i = 1:surf.ndiv-1
                #     delu[i] = deltacurve(suc[i])
                # end
                # println("Remeshed delta")
                
                delu[:] = del_iter[:]
                Eu[:] = E_iter[:]
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
                #srcstr = (del_iter[end] - del_iter[end-1])/sqrt(Re)
                #if srcstr > 0.
                #    srcstr = 0.
                #end
                #println(srcstr)
                #push!(wk_src, TwoDSource(surf.bnd_x_u[surf.ndiv] + 0.5*surf.uref*dt, 0.5*(surf.bnd_z_u[surf.ndiv] + surf.bnd_z_l[surf.ndiv]), srcstr*20))
            end

            if iter > iterMax
                println("failed to converge")
                delu[:] = del_iter[:]
                Eu[:] = E_iter[:]
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
                break
            end
            
            if iter == 3
                #plot(surf.x, qu)
                plot(surf.x[2:end], del_iter)
                #println(del_iter)
                #  error("here")
            end
            
        end        
        
        #Calculate adot
        update_atermdot(surf, dt)
        
        #Set previous values of aterm to be used for derivatives in next time step
        surf.a0prev[1] = surf.a0[1]
        for ia = 1:3
            surf.aprev[ia] = surf.aterm[ia]
        end

        #Update induced velocities to include effect of last shed vortex
        update_indbound(surf, curfield)
        
        #Calculate bound vortex strengths
        update_bv_src(surf)
        
        #Wake rollup
        wakeroll(surf, curfield, dt)
        
        #Force calculation
        cnc, cnnc, cn, cs, cl, cd, int_wax, int_c, int_t = calc_forces(surf, int_wax, int_c, int_t, dt)
        
        qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl)
            end
        end
        
        #LE velocity and stagnation point location
        
        vle = qu[1]
        
        if vle > 0.
            qspl = Spline1D(surf.x, ql)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        else
            qspl = Spline1D(surf.x, qu)
            stag = try
                roots(qspl, maxn=1)[1]
            catch
                0.
            end
        end
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, vle,
                        cl, cd, cnc, cnnc, cn, cs, stag])

    end    

    mat = mat'

    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, mat)
    close(f)

    return mat, surf, curfield

end


mat, surf, curfield = IBLnew(surf, curfield, nsteps, dtstar, startflag, writeflag, writeInterval, delvort)


