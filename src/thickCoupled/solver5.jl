function IBL_shape_attached(Re, surf::TwoDSurfThick, curfield::TwoDFlowField, nsteps::Int64 = 300, dtstar::Float64 = 0.015, startflag = 0, writeflag = 0, writeInterval = 1000., delvort = delNone(); maxwrite = 50, nround=6)
    
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

    #At this time coded for symmetric situations only
    qu = zeros(surf.ndiv)
    
    s = zeros(surf.ndiv)
    quc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)
    qu_prev = zeros(surf.ndiv)
    qucx = zeros(surf.ndiv-1)
    quct = zeros(surf.ndiv-1)
    sc = zeros(surf.ndiv-1)

    ql =  zeros(surf.ndiv) 
    qlc = zeros(surf.ndiv-1)
    qlc_prev = zeros(surf.ndiv-1)
    ql_prev = zeros(surf.ndiv)
    qlcx = zeros(surf.ndiv-1)
    qlct = zeros(surf.ndiv-1)

    thick_orig = zeros(surf.ndiv)
    thick_orig[:] = surf.thick[:]
    thick_slope_orig = zeros(surf.ndiv)
    thick_slope_orig[:] = surf.thick_slope[:]

    cam_orig = zeros(surf.ndiv)
    cam_slope_orig = zeros(surf.ndiv)
    cam_orig[:] = surf.cam[:]
    cam_slope_orig[:] = surf.cam_slope[:]

   #Initialise boundary layer
    delu,dell, Eu, El = initDelE(surf.ndiv-1)

    vle = 0.
    stindex = 0.
    
    dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    s[1] = 0.
    for i = 2:surf.ndiv
        s[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
    end

    #surf.su_i[:] =thick_slope_orig[:]


    for i = 1:surf.ndiv-1
        sc[i] = (s[i] + s[i+1])/2
    end

    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)
    cpu =   zeros(surf.ndiv)
    cpl =  zeros(surf.ndiv)

    #wk_src = TwoDSource[]

    x_w = collect(surf.c*1.01:surf.c*0.01:surf.c*3)
    nw = length(x_w)
    wfn = zeros(nw)
    
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
        delu_iter = zeros(surf.ndiv-1)
        delu_prev = zeros(surf.ndiv-1)
   
	dell_iter = zeros(surf.ndiv-1)
	dell_prev = zeros(surf.ndiv-1)

	El_iter = zeros(surf.ndiv-1)
	Eu_iter = zeros(surf.ndiv-1)

        soln = zeros(2*surf.naterm+1)

    	#stindex = 0.
        #stindex_prev = 0.
        resTol = 1e-5
        iterMax = 10

        quc_prev[:] = quc[:]       
        qlc_prev[:] = qlc[:]       
#	qu_prev[:] = qu[:]
#	ql_prev[:] = ql[:]
        
        #while (iter < iterMax)
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

	    qucavg = (quc+qlc)/2
	    delavg = (delu_iter + dell_iter)/2

            a_wfn = log10(2*qucavg[end]*delavg[end]/sqrt(Re))
            wfn[1] = 2*(qucavg[end]*delavg[end] - qucavg[end-1]*delavg[end-1])/(sqrt(Re)*(surf.x[end] - surf.x[end-1]))
            for i = 2:nw
             #  wfn[i] = (10^(a_wfn - 3.2*(x_w[i] - surf.c)) - 10^(a_wfn - 3.2*(x_w[i] - 1.)))/(x_w[i] - x_w[i-1])
            	wfn[i] = 10^(a_wfn - 10^(a_wfn - 3.2*(x_w[i] - 1.)))
            end
            
            #Source strengths in wake and induced velocity
            uind_src = zeros(surf.ndiv)
            for i = 1:surf.ndiv
                for iw = 1:nw-1
                    str = 0.5(wfn[iw] + wfn[iw+1])
                    xloc = 0.5*(x_w[iw] + x_w[iw+1])
                    uind_src[i] += 1/(2*pi)*str/(xloc - surf.x[i])*(x_w[iw+1] - x_w[iw])
                end
            end

            surf.uind_u[:] += uind_src[:]
            surf.uind_l[:] += uind_src[:]
            
            ind_new_u_u, ind_new_w_u = ind_vel([temptev], surf.bnd_x_u, surf.bnd_z_u)
            ind_new_u_l, ind_new_w_l = ind_vel([temptev], surf.bnd_x_l, surf.bnd_z_l)
            
            surf.uind_u[:] += ind_new_u_u[:]
            surf.wind_u[:] += ind_new_w_u[:]
            surf.uind_l[:] += ind_new_u_l[:]
            surf.wind_l[:] += ind_new_w_l[:]

	     qu[:], ql[:], phi_u[:], phi_l[:], cpu[:], cpl[:] = calc_edgeVel_cp(surf, [curfield.u[1], curfield.w[1]], phi_u, phi_l, dt)

     
            surf.uind_u[:] -= uind_src[:]
            surf.uind_l[:] -= uind_src[:]
          
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]
            
            smoothScaledEnd!(surf.x, qu,10)
            smoothScaledEnd!(surf.x, ql,10)           

	vle = qu[1]
	#surf.su[:] = su[:]

	surf.ueU[:] = qu[:] 
	surf.ueL[:] = ql[:] 

	if vle > 0.
		
	     stindex = argmin(ql)
	    # println("stagnation index ",stindex)
	     #println(stindex)
	     if stindex == 140
		    error("stop")
	     end
	     
	     x_u, x_l, su, sl, qustag, qlstag, qustag_prev, qlstag_prev, delustag, dellstag, Eustag, Elstag, dellstag_iter, delustag_iter, delustag_prev, dellstag_prev, Eustag_iter, Elstag_iter = reconstructGrid(stindex, surf, s, qu, ql, quc_prev, qlc_prev, delu, dell, Eu, El, delu_iter, dell_iter, delu_prev, dell_iter, Eu_iter, El_iter)

	     upper_size = length(x_u) 
	     lower_size = length(x_l)
	     #println("maximum upper ",maximum(x_u))

	     qustagc = zeros(upper_size-1)
	     qlstagc = zeros(lower_size-1)

	     qustagc_prev = zeros(upper_size-1)	
	     qlstagc_prev = zeros(lower_size-1)


	     qustagx =  zeros(upper_size-1)
	     qlstagx = zeros(lower_size-1)

	     qustagct =  zeros(upper_size-1)
	     qlstagct = zeros(lower_size-1)

	     suc =  zeros(upper_size-1)
	     slc = zeros(lower_size-1)

	     suc[1:end] = (su[2:end] + su[1:end-1])./2
	     slc[1:end] = (sl[2:end] + sl[1:end-1])./2

	     qustagc_prev[:] = qustag_prev[:] 	
	     qlstagc_prev[:] = qlstag_prev[:]	

	 else
		stindex = argmin(qu) 
		qustag = qu[stindex+1:end] 
		qlstag = [reverse(qu[1:stindex]);ql[stindex+1:end]]
	       
		qustlen = length(qustag)
		qlstlen = length(qlstag)
		    
		surf.qlstlen = qlstlen
		surf.qustlen = qustlen

		xustag = zeros(qustlen)
		xlstag = zeros(qlstlen)
		error("upper stagnation point on a postive aoa")

	 end


            #Solve the FV problem at cell centres

           # for i = 1:surf.ndiv-1
             quc[1:end] = (qu[2:end] + qu[1:end-1])/2
	     qlc[1:end] = (ql[2:end] + ql[1:end-1])/2
	     qustagc[1:end] = (qustag[2:end] + qustag[1:end-1])./2
	     qlstagc[1:end] = (qlstag[2:end] + qlstag[1:end-1])./2
           # end
            

	    qustagx[1:end] = diff1(suc, qustagc)
	   #qustagx[2:end] = diff(qustagc)./diff(suc)
           #qustagx[1] = qustagx[2]

	    qlstagx[1:end] = diff1(slc, qlstagc)
	   #qlstagx[2:end] = diff(qlstagc)./diff(slc)
           #qlstagx[1] = qlstagx[2]

            #smoothEnd!(qucx, 10)
           
	   # surf.ueU[:] = qu[:] 
	   # surf.ueL[:] = ql[:] 

            if istep == 1
		    qustagct[:] .= 0.
		    qlstagct[:] .= 0.
            else
             
		     qustagct[:] = (qustagc[:] - qustagc_prev[:])/dt
		     qlstagct[:] = (qlstagc[:] - qlstagc_prev[:])/dt

            end

            wu = [delustag delustag.*(Eustag .+ 1.0)]
            wl = [dellstag dellstag.*(Elstag .+ 1.0)]

            wusoln, i_sepu = FVMIBLgridvar(wu, qustagc, qustagct, qustagx, diff(su), t-dt, t)
            wlsoln, i_sepl = FVMIBLgridvar(wl, qlstagc, qlstagct, qlstagx, diff(sl), t-dt, t)
	    
	    delustag_prev[:] = delustag_iter[:]
            delustag_iter[:] = wusoln[:,1]

            dellstag_prev[:] = dellstag_iter[:]
            dellstag_iter[:] = wlsoln[:,1]

	    Eustag_iter[:] = wusoln[:,2]./wusoln[:,1] .- 1.
            Elstag_iter[:] = wlsoln[:,2]./wlsoln[:,1] .- 1

            #smoothScaledEnd!(suc, delustag_iter, 10)
            #smoothScaledEnd!(slc, dellstag_iter, 10)	    
	    

	 quc, qlc, quc_prev, qlc_prev, delu, dell, Eu, El, delu_iter, dell_iter, delu_prev, dell_iter, Eu_iter, El_iter = reverseReconstructGrid(stindex, surf, qustag, qlstag,  qustagc, qlstagc, qustagc_prev, qlstagc_prev, delustag, dellstag, Eustag, Elstag, dellstag_iter, delustag_iter, delustag_prev, dellstag_prev, Eustag_iter, Elstag_iter)
	   
	   #figure("dell")
	   #plot(dellstag_iter)		
	    #Find suitable naca coefficients to fit the modified airfoil
	    smoothScaledEnd!(sc, delu_iter, 10)
	    smoothScaledEnd!(sc, dell_iter, 10)	    

            newthick = zeros(surf.ndiv)
 	    conU = zeros(surf.ndiv)
	    conL = zeros(surf.ndiv)
	    thickUpdate = zeros(surf.ndiv)
	    
	    newcamb = zeros(surf.ndiv)
	    cambConU = zeros(surf.ndiv)
	    cambConL = zeros(surf.ndiv)
	    camUpdate = zeros(surf.ndiv)

            
	    for i = 1:surf.ndiv-1

	    	conU[i] =  (quc[i]*delu[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i] + cam_slope_orig[i]).^2))
	    	conL[i] =  (qlc[i]*dell[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i] + cam_slope_orig[i]).^2))
                
		thickUpdate[i] = (conU[i] + conL[i])/2
		camUpdate[i] = (conU[i] - conL[i])/2.0
		
		newcamb[i] = cam_orig[i] +  camUpdate[i]
		newthick[i] = thick_orig[i] + thickUpdate[i] 

            end
            
	    conU[surf.ndiv] =  (quc[surf.ndiv-1]*delu[surf.ndiv-1])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[surf.ndiv] + cam_slope_orig[surf.ndiv]).^2)) 
	    conL[surf.ndiv] = (qlc[surf.ndiv-1]*dell[surf.ndiv-1])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[surf.ndiv] + cam_slope_orig[surf.ndiv]).^2))  

            newthick[surf.ndiv] = thick_orig[surf.ndiv] + 0.5*(conU[surf.ndiv] + conL[surf.ndiv]) 
	    newcamb[surf.ndiv] = cam_orig[surf.ndiv] + 0.5*(conU[surf.ndiv] - conL[surf.ndiv]) 

	    surf.su_i[:] = newcamb[:] 

	    bstart = [-0.1260; -0.3516; 0.2843; -0.1015; 0.0; 0.0]
	    bstart_cam = [0.0; 0.0; 0.0; 0.0] 

            coef = find_nacaCoef(surf, newthick, bstart)
	    coef_cam = find_nacaCamb(surf, newcamb, bstart_cam)

            th = parse(Int, surf.coord_file[7:8])/100.
            b1 = 0.2969
            b = [b1; coef]
	    c = coef_cam
            
	    @. nacath(x) = 5*th*(b[1]*sqrt(x) + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4 + b[6]*x^5 + b[7]*x^6)
	    @. nacacam(x) = (c[1] + c[2]*x + c[3]*x^2 + c[4]*x^3)
            
            #Find new shape of airfoil
            for i = 1:surf.ndiv
                surf.thick[i] = nacath(surf.x[i])
                surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
		
		surf.cam[i] = nacacam(surf.x[i])
		surf.cam_slope[i] = ForwardDiff.derivative(nacacam, surf.x[i])
            end
	    println(iter, "   ", res," time: ",t, " stindex: ", stindex)
	   # if iter ==2
	   # println(c)
	   #error("stop here, checkpoint")
           # end
            #Check for convergence
            res =  sum(abs.(delu_prev .- delu_iter))/length(delu_prev) #+ sum(abs.(dellstag_prev .- dellstag_iter))sum(abs.(resU))/norm(x_u) 
            resL =   sum(abs.(dell_prev[10:end] .- dell_iter[10:end]))

            #if iter == iterMax
            if res <= resTol
                println("converged")
		#println("Lower res ", resL)
                delu[:] = delu_iter[:]
                Eu[:] = Eu_iter[:]

                dell[:] = dell_iter[:]
                El[:] = El_iter[:]
		

		surf.deltaU[2:end] = delu[:]
		#surf.ueU[:] = qu[:] 
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            end

            if iter == 3 && mod(istep,10) == 0
                figure("Edge velocity")
		clf()
                plot(surf.x, qu) 
		plot(surf.x, ql)
                figure("Thickness")
                plot(surf.x, surf.thick)
                axis("equal")


		clf()
                figure("delta Higher")
                plot(suc, delustag_iter)
		figure("delta Lower")
                plot(slc, dellstag_iter)
		
		#figure("E Higher")
		#plot(surf.x[2:end], Eu)
	#	figure("E Lower")
		#plot(surf.x[2:end], El)

		#error("first plot")
            end

	    if istep ==25
		@bp	
		#error("stop")
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
        
        qu[:], ql[:], phi_u[:], phi_l[:], cpu[:], cpl[:] = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)
        
        # write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(t,sigdigits=nround))"
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl, suc, delu, Eu, thick_orig, quc, qucx, quct)
            end
        end
        
        #LE velocity and stagnation point location
        
        
        mat = hcat(mat,[t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u,
                        cl, cd, cnc, cnnc, cn, cs, vle, stindex])
    end    

    mat = mat'

    f = open("resultsSummary", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, mat)
    close(f)

    return mat, surf, curfield

end

function diff1(x::Array{Float64,1}, f::Array{Float64,1})

    fp = zeros(length(x))
    fpp = zeros(length(x))
    dx = x[2:end] - x[1:end-1]
    df = f[2:end] - f[1:end-1]

    fpp[1:end-1] = atan.(df,dx)

    dx1 = x[2:end-1] - x[1:end-2]
    dx2 = x[3:end] - x[2:end-1]
    ang = (dx2.*fpp[1:end-2] + dx1.*fpp[2:end-1])./(dx1 .+ dx2)
    fp[2:end-1] = tan.(ang)

    fp[1] = 2.0*tan.(fpp[1])- fp[2]
    fp[end] = 2.0*tan.(fpp[end-1])- fp[end-1]

    return fp


end

