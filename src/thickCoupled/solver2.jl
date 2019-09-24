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
    
    su = zeros(surf.ndiv)
    quc = zeros(surf.ndiv-1)
    quc_prev = zeros(surf.ndiv-1)
    qu_prev = zeros(surf.ndiv)
    qucx = zeros(surf.ndiv-1)
    quct = zeros(surf.ndiv-1)
    suc = zeros(surf.ndiv-1)

    ql =  zeros(surf.ndiv) 
    #qlc = zeros(surf.ndiv-1)
    qluc_prev = zeros(surf.ndiv-1)
    ql_prev = zeros(surf.ndiv)
    qlcx = zeros(surf.ndiv-1)
    qlct = zeros(surf.ndiv-1)

    thick_orig = zeros(surf.ndiv)
    thick_orig[:] = surf.thick[:]
    thick_slope_orig = zeros(surf.ndiv)
    thick_slope_orig[:] = surf.thick_slope[:]

   #Initialise boundary layer
    #delu,dell, Eu, El = initDelE(surf.ndiv-1)

    
    dsdx = zeros(surf.ndiv)
    for i = 2:surf.ndiv
        dsdx[i] = sqrt(1 + (surf.cam_slope[i] + surf.thick_slope[i])^2)
    end
    dsdx[1] = dsdx[2]
    su[1] = 0.
    for i = 2:surf.ndiv
        su[i] = simpleTrapz(dsdx[1:i], surf.x[1:i])
    end

    surf.su_i[:] =thick_slope_orig[:]


    for i = 1:surf.ndiv-1
        suc[i] = (su[i] + su[i+1])/2
    end

    phi_u = zeros(surf.ndiv)
    phi_l = zeros(surf.ndiv)

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
        del_iter = zeros(surf.ndiv-1)
        del_prev = zeros(surf.ndiv-1)
        soln = zeros(2*surf.naterm+1)
        E_iter = zeros(surf.ndiv-1)

        resTol = 1e-5
        iterMax = 10

        quc_prev[:] = quc[:]       
        qlc_prev[:] = qlc[:]       
	qu_prev[:] = qu[:]
	ql_prev[:] = ql[:]
        
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

            a_wfn = log10(2*quc[end]*del_iter[end]/sqrt(Re))
            wfn[1] = 2*(quc[end]*del_iter[end] - quc[end-1]*del_iter[end-1])/(sqrt(Re)*(surf.x[end] - surf.x[end-1]))
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

            #qu[:], ql[:] = calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])
qu, ql, phi_u, phi_l, cpu, cpl = calc_edgeVel_cp(surf, [curfield.u[1]; curfield.w[1]], phi_u, phi_l, dt)

     
            vle = qu[1]
	    surf.su[:] = su[:]

	    surf.ueU[:] = qu[:] 
	    surf.ueL[:] = ql[:] 
        
    	    if vle > 0.
        	    
		   stindex = argmin(ql)
		   println("stagnation index ",stindex)
		   println(stindex) 
		   x_u, x_l, su, sl, qustag, qlstag, qustag_prev, qlstag_prev  = reconstructGrid(stindex, surf, su, qu, ql, qu_prev, ql_prev)

		
		   upper_size = length(x_u) 
		   lower_size = length(x_l)
		   println("maximum upper ",maximum(x_u))
	#	   println("upper size ",upper_size)
		   
		   qustagc = zeros(upper_size-1)
		   qlstagc = zeros(lower_size-1)
		   
		   qucstag_prev = zeros(upper_size-1)	
		   qlcstag_prev = zeros(lower_size-1)


		   qustagx =  zeros(upper_size-1)
		   qlstagx = zeros(lower_size-1)

		   qustagct =  zeros(upper_size-1)
		   qlstagct = zeros(lower_size-1)

                   delu_iter = zeros(upper_size-1)
                   delu_prev = zeros(upper_size-1)

                   dell_iter = zeros(lower_size-1)
                   dell_prev = zeros(lower_size-1)
		   
		   El_iter = zeros(lower_size-1)
	           Eu_iter = zeros(upper_size-1)


		   suc =  zeros(upper_size-1)
		   slc = zeros(lower_size-1)
		  
		   suc[1:end] = (su[2:end] + su[1:end-1])./2
		   slc[1:end] = (sl[2:end] + sl[1:end-1])./2

		   qustagc_prev[1:end] = (qustag_prev[1:end-1] + qustag_prev[2:end])./2	
		   qlstagc_prev[1:end] = (qlstag_prev[1:end-1] + qlstag_prev[2:end])./2	
		   delu,dell, Eu, El = initDelE(upper_size-1, lower_size-1)	
       	     
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

	     end

	  # qustag = [reverse(ql[1:stindex]);qu] 
	 #  qlstag = ql[stindex+1:end]

            surf.uind_u[:] -= uind_src[:]
            surf.uind_l[:] -= uind_src[:]
          
            surf.uind_u[:] -= ind_new_u_u[:]
            surf.wind_u[:] -= ind_new_w_u[:]
            surf.uind_l[:] -= ind_new_u_l[:]
            surf.wind_l[:] -= ind_new_w_l[:]
            
            smoothScaledEnd!(x_u, qustag,10)
            smoothScaledEnd!(x_l, qlstag,10)           

            #Solve the FV problem at cell centres

           # for i = 1:surf.ndiv-1
             quc[1:end] = (qu[2:end] + qu[1:end-1])/2
	     qlc[1:end] = (ql[2:end] + ql[1:end-1])/2
	     qustagc[1:end] = (qustag[2:end] + qustag[1:end-1])./2
	     qlstagc[1:end] = (qlstag[2:end] + qlstag[1:end-1])./2
           # end
            
           
	   qustagx[2:end] = diff(qustagc)./diff(suc)
           qustagx[1] = qustagx[2]

	   qlstagx[2:end] = diff(qlstagc)./diff(slc)
           qlstagx[1] = qlstagx[2]

            #smoothEnd!(qucx, 10)
           
	   # surf.ueU[:] = qu[:] 
	   # surf.ueL[:] = ql[:] 

	   figure("velocity")
	   plot(x_u, qustag, label="upper")
	   #figure("lower stagnation")
	   plot(x_l, qlstag, label="lower")
	  # legend()
	  # surf.su[:] = su[:]

	   error("stop here")
            if istep == 1
             qustagct[:] .= 0.
	     qlstagct[:] .= 0.
            else
             
	     qustagct[:] = (qustagc[:] - qustagc_prev[:])/dt
 	     qlstagct[:] = (qlstagc[:] - qlstagc_prev[:])/dt

            end

            wu = [delu delu.*(Eu .+ 1)]
            wl = [dell dell.*(El .+ 1)]

           wusoln, i_sepu = FVMIBLgridvar(wu, qustagc, qustagct, qustagx, diff(su), t-dt, t)
           wlsoln, i_sepl = FVMIBLgridvar(wl, qlstagc, qlstagct, qlstagx, diff(sl), t-dt, t)

            delu_prev[:] = delu_iter[:]
            delu_iter[:] = wusoln[:,1]

            dell_prev[:] = dell_iter[:]
            dell_iter[:] = wlsoln[:,1]

	    Eu_iter[:] = wusoln[:,2]./wusoln[:,1] .- 1.
            El_iter[:] = wlsoln[:,2]./wlsoln[:,1] .- 1

            smoothScaledEnd!(suc, delu_iter, 10)
            smoothScaledEnd!(slc, dell_iter, 10)
	    
	    deluu, delll, quc, qlc = reverseReconstructGrid(stindex, delu_iter, dell_iter,qustagc, qlstagc)   


            
	    #Find suitable naca coefficients to fit the modified airfoil

            newthick = zeros(surf.ndiv)
 	    thickconU = zeros(surf.ndiv)
	    thickconL = zeros(surf.ndiv)
	    thickUpdate = zeros(surf.ndiv)

            
	    for i = 1:surf.ndiv-1

	    	thickconU[i] =  (quc[i]*deluu[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i]).^2))
	    	thickconL[i] =  (qlc[i]*delll[i])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[i]).^2))
                thickUpdate[i] = (thickconU[i] + thickconL[i])/2
		newthick[i] = thick_orig[i] + thicjUpdate[i] 
            end
            
            newthick[surf.ndiv] = thick_orig[surf.ndiv] + (quc[surf.ndiv-1]*del_iter[surf.ndiv-1])/(sqrt(Re)*sqrt(1. + (thick_slope_orig[surf.ndiv]).^2))

            surf.thick_a[:] = newthick[:]
            bstart = [-0.1260; -0.3516; 0.2843; -0.1015; 0.0; 0.0; 0.0; 0.0]

            coef = find_nacaCoef(surf, newthick, bstart)

            th = parse(Int, surf.coord_file[7:8])/100.
            b1 = 0.2969
            b = [b1; coef]
            
            @. nacath(x) = 5*th*(b[1]*sqrt(x) + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4 + b[6]*x^5 + b[7]*x^6 + b[8]*x^7 + b[9]*x^8)
            
            #Find new shape of airfoil
            for i = 1:surf.ndiv
                surf.thick[i] = nacath(surf.x[i])
                surf.thick_slope[i] = ForwardDiff.derivative(nacath, surf.x[i])
            end
            
            #Check for convergence
            res =  sum(abs.(del_prev .- del_iter))
            println(iter, "   ", res," time: ",t)


            #if iter == iterMax
            if res <= resTol
                println("converged")
                delu[:] = del_iter[:]
                Eu[:] = E_iter[:]
		surf.deltaU[2:end] = delu[:]
		#surf.ueU[:] = qu[:] 
                push!(curfield.tev, TwoDVort(xloc_tev, zloc_tev, tevstr, vcore, 0., 0.))
            end

            if iter == 3 && mod(istep,10) == 0
                figure(1)
                plot(surf.x, qu)
                figure(2)
                plot(surf.x, surf.thick)
                axis("equal")
                figure(3)
                plot(surf.x[2:end], delu)
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
                writeStamp(dirname, t, surf, curfield, qu, ql, cpu, cpl, suc, delu, Eu, thick_orig, quc, qucx, quct)
            end
        end
        
        #LE velocity and stagnation point location
        
        
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


