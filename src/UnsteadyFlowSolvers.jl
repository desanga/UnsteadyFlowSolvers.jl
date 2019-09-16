#__precompile__(true)

module UnsteadyFlowSolvers

import Dierckx: Spline1D, derivative, evaluate, roots

using ForwardDiff

import DelimitedFiles: writedlm, readdlm

import Serialization

import NLsolve: nlsolve, not_in_place

import Statistics: mean

import PyPlot: plot, scatter, figure, xlabel, ylabel, xlim, ylim,
    xticks, yticks, subplot, subplot2grid, legend, axis, savefig,
    close, tight_layout, clf

import SpecialFunctions: besselh, besselk

import Plots: @layout

import LaTeXStrings: @L_str

#For use in development and debugging
import Revise

using Printf

using DataStructures

using Interpolations

import QuadGK: quadgk

import LinearAlgebra: det

import Roots: find_zero, Bisection

using LeastSquaresOptim

import Optim: optimize, BFGS

import LsqFit: curve_fit

using LsqFit

export
    # kinematics types and funtions
    MotionDef,
    KinemPar,
    KinemDef,
    EldUpDef,
    EldDownDef,
    EldUpDownDef,
    EldUptstartDef,
    ConstDef,
    EldRampReturnDef,
    EldUpIntDef,
    EldUpInttstartDef,
    SinDef,
    CosDef,

    # 2D low-order solver types
    TheoDef,
    TwoDSurf,
    TwoDOFPar,
    KinemPar2DOF,
    TwoDVort,
    TwoDFlowField,
    KelvinCondition,
    KelvinKutta,
    TwoDSurfThick,

    # vortex count control utility
    delVortDef,
    delNone,
    delSpalart,

    # utility functions
    simpleTrapz,
    camber_calc,
    find_tstep,
    simpleInterp,
    cleanWrite,
    getEndCycle,
    secant_method,
    errorCalc,
    smoothEdges!,
    smoothEnd!,
    smoothStart!,

    #cfd utility function
    compareForces,
    compare_cp_edgevel,

    # 2D low-order solver methods
    theodorsen,
    lautat,
    ldvm,
    ldvmLin,
    ldvm2DOF,

    # 2D plot output functions
    makeForcePlots2D,
    makeVortPlots2D,

    # 2D postprocessing functions
    calc_edgeVel,
    calc_edgeVel_cp,

    #IBLThickCoupled,
    initDelE,
    FVMIBLorig,
    dtfunFVM,
    dtfunjac,
    transpTogetherWake,
    transpSimul,
    transResidual,
    FVMIBLgridvar,
    IBLshed,
    IBLsheddisp,
    theodorsen_stag,
    qs_utat,
    theodorsen_delcp,
    theodorsen_stag_nonlinear,
    calc_stag,
    IBL_simul,
    IBL_simul_iter,
    IBL_simul_iter_shed,
    update_thickLHS,
    update_thickRHS,
    update_kinem,
    update_externalvel,
    update_boundpos,
    update_indbound,
    ind_vel,
    update_atermdot,
    update_bv_src,
    wakeroll,
    calc_forces,
    writeStamp,
    smoothEnd!,
    TwoDSource,
    ind_vel_src,
    IBL_shape_attached,
    IBL_vel_attached,
    find_nacaCoef,
    smoothScaledEnd!



### source files

# kinematic types
include("kinem.jl")

# utility functions
include("utils.jl")
include("cfdutils.jl")
# vortex count control utility
include("delVort.jl")

# low-order 2D thin solvers
include("thin2D/typedefs.jl")            # type definitions
include("thin2D/calcs.jl")               # calculation functions
include("thin2D/solvers.jl")             # solver methods
include("thin2D/postprocess.jl")         # postprocessing functions

# low-order 2D thick solvers
include("thick2D/typedefs.jl")            # type definitions
include("thick2D/calcs.jl")               # calculation functions
include("thick2D/solvers.jl")             # solver methods
include("thick2D/postprocess.jl")         # postprocessing functions


# low-order 2D thick IBL solvers
#include("thickCoupled/typedefs.jl")            # type definitions
include("thickCoupled/calcs.jl")               # calculation functions
include("thickCoupled/solver1.jl")             # solver methods
include("thickCoupled/postprocess.jl")         # postprocessing functions

# 2D plotting functions
include("plots/plots2D.jl")

# IBL thick-coupled functions
#include("thickCoupled/IBLThickCoupled.jl")
#include("thickCoupled/IBLFV.jl")
#include("thickCoupled/correlate.jl")
#include("thickCoupled/transpCoupled.jl")

# IBL thin-coupled functions
#include("thinCoupled/IBLCoupled.jl")
#include("thinCoupled/IBLFV.jl")
#include("thinCoupled/correlate.jl")


end
