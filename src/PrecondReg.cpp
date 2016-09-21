#ifndef _PRECONDREG_CPP_
#define _PRECONDREG_CPP_

#include "PrecondReg.hpp"



namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PrecondReg"
PrecondReg::PrecondReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~PrecondReg"
PrecondReg::~PrecondReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "OptimalControlRegistration"
PrecondReg::PrecondReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode PrecondReg::Initialize()
{
    PetscFunctionBegin;

    this->m_Opt = NULL; ///< options (default; on grid we solve)
    this->m_OptCoarse = NULL; ///< options for coarse grid

    this->m_MatVec = NULL; ///< pointer to matvec in krylov method
    this->m_MatVecEigEst = NULL; ///< pointer to matvec in krylov method
    this->m_KrylovMethod = NULL; ///< pointer to krylov method
    this->m_KrylovMethodEigEst = NULL; ///< pointer to krylov method
    this->m_RandomNumGen = NULL; ///< pointer to krylov method

    this->m_xCoarse = NULL; ///< container for input to hessian matvec on coarse grid
    this->m_HxCoarse = NULL; ///< container for hessian matvec on coarse grid

    this->m_PreProc = NULL; ///< pointer to preprocessing operator
    this->m_OptProbCoarse = NULL; ///< optimization problem on coarse grid

    this->m_ControlVariable = NULL; ///< control variable on fine grid
    this->m_IncControlVariable = NULL; ///< incremental control variable on fine grid

    this->m_StateVariableCoarse = NULL; ///< state variable on coarse grid
    this->m_AdjointVariableCoarse = NULL; ///< adjoint variable on coarse grid
    this->m_ControlVariableCoarse = NULL; ///< control variable on coarse grid
    this->m_IncControlVariableCoarse = NULL; ///< incremental control variable on coarse grid

    this->m_WorkVecField = NULL; ///< temporary vector field
    this->m_WorkScaField1 = NULL; ///< temporary scalar field
    this->m_WorkScaField2 = NULL; ///< temporary scalar field
    this->m_WorkScaFieldCoarse1 = NULL; ///< temporary scalar field (coarse level)
    this->m_WorkScaFieldCoarse2 = NULL; ///< temporary scalar field (coarse level)

    this->m_EigenValuesEstimated=false;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PrecondReg::ClearMemory()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_KrylovMethod != NULL) {
        ierr = KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }
    if (this->m_MatVec != NULL) {
        ierr = MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }
    if (this->m_MatVecEigEst != NULL) {
        ierr = MatDestroy(&this->m_MatVecEigEst); CHKERRQ(ierr);
        this->m_MatVecEigEst = NULL;
    }


    if (this->m_xCoarse != NULL) {
        ierr = VecDestroy(&this->m_xCoarse); CHKERRQ(ierr);
        this->m_xCoarse = NULL;
    }
    if (this->m_HxCoarse != NULL) {
        ierr = VecDestroy(&this->m_HxCoarse); CHKERRQ(ierr);
        this->m_HxCoarse = NULL;
    }
    if (this->m_StateVariableCoarse != NULL) {
        ierr = VecDestroy(&this->m_StateVariableCoarse); CHKERRQ(ierr);
        this->m_StateVariableCoarse = NULL;
    }
    if (this->m_AdjointVariableCoarse != NULL) {
        ierr = VecDestroy(&this->m_AdjointVariableCoarse); CHKERRQ(ierr);
        this->m_AdjointVariableCoarse = NULL;
    }


    if (this->m_WorkVecField != NULL) {
        delete this->m_WorkVecField;
        this->m_WorkVecField=NULL;
    }
    if (this->m_WorkScaField1 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField1); CHKERRQ(ierr);
        this->m_WorkScaField1 = NULL;
    }
    if (this->m_WorkScaField2 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaField2); CHKERRQ(ierr);
        this->m_WorkScaField2 = NULL;
    }
    if (this->m_WorkScaFieldCoarse1 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaFieldCoarse1); CHKERRQ(ierr);
        this->m_WorkScaFieldCoarse1 = NULL;
    }
    if (this->m_WorkScaFieldCoarse2 != NULL) {
        ierr = VecDestroy(&this->m_WorkScaFieldCoarse2); CHKERRQ(ierr);
        this->m_WorkScaFieldCoarse2 = NULL;
    }


    if (this->m_OptProbCoarse != NULL) {
        delete this->m_OptProbCoarse;
        this->m_OptProbCoarse=NULL;
    }

    if (this->m_ControlVariable != NULL) {
        delete this->m_ControlVariable;
        this->m_ControlVariable=NULL;
    }
    if (this->m_IncControlVariable != NULL) {
        delete this->m_IncControlVariable;
        this->m_IncControlVariable=NULL;
    }


    if (this->m_ControlVariableCoarse != NULL) {
        delete this->m_ControlVariableCoarse;
        this->m_ControlVariableCoarse=NULL;
    }
    if (this->m_IncControlVariableCoarse != NULL) {
        delete this->m_IncControlVariableCoarse;
        this->m_IncControlVariableCoarse=NULL;
    }

    if (this->m_OptCoarse != NULL) {
        delete this->m_OptCoarse;
        this->m_OptCoarse = NULL;
    }

    if (this->m_RandomNumGen!=NULL) {
        ierr = PetscRandomDestroy(&this->m_RandomNumGen); CHKERRQ(ierr);
        this->m_RandomNumGen=NULL;
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetProblem"
PetscErrorCode PrecondReg::SetProblem(PrecondReg::OptProbType* optprob) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_OptProb = optprob;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetPreProc"
PetscErrorCode PrecondReg::SetPreProc(PreProcReg* preproc) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr = Assert(preproc!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_PreProc = preproc;

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief make sure that we re-initiate/recompute important
 * quantities; implemented to allow multiple calls of the solver
 * without destroying it;
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Reset"
PetscErrorCode PrecondReg::Reset() {
    PetscErrorCode ierr = 0;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype) {
        case NOPC:
        {
            // no need to do anything
            break;
        }
        case INVREG:
        {
            // no need to do anything
            break;
        }
        case TWOLEVEL:
        {
            // in case we call the solver multiple times (for
            // instance when we do parameter continuation) without
            // destroying the preconditioner, it is necessary to
            // recompute eigenvalues
            this->m_EigenValuesEstimated = false;
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief setup phase of preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DoSetup"
PetscErrorCode PrecondReg::DoSetup() {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // start timer
    ierr = this->m_Opt->StartTimer(PMVSETUP); CHKERRQ(ierr);

    // switch case for choice of preconditioner
    if (this->m_Opt->GetKrylovSolverPara().pctype == TWOLEVEL) {
        ierr = this->Setup2LevelPrecond(); CHKERRQ(ierr);
    }
    this->m_Opt->PrecondSetupDone(true);

    // stop timer
    ierr = this->m_Opt->StopTimer(PMVSETUP); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MatVec"
PetscErrorCode PrecondReg::MatVec(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype) {
        case NOPC:
        {
            ierr = WrngMsg("no preconditioner used"); CHKERRQ(ierr);
            ierr = VecCopy(x,Px); CHKERRQ(ierr);
            break;
        }
        case INVREG:
        {
            ierr = this->ApplyInvRegPC(Px, x); CHKERRQ(ierr);
            break;
        }
        case TWOLEVEL:
        {
            ierr = this->Apply2LevelPC(Px, x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief apply inverse of regularization operator as preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegPC"
PetscErrorCode PrecondReg::ApplyInvRegPC(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // check if optimization problem is set up
    ierr = Assert(this->m_OptProb != NULL, "null pointer"); CHKERRQ(ierr);

    // start timer
    ierr = this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // apply inverse regularization operator
    ierr = this->m_OptProb->ApplyInvRegOp(Px, x); CHKERRQ(ierr);

    // stop timer
    ierr = this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Apply2LevelPC"
PetscErrorCode PrecondReg::Apply2LevelPC(Vec Px, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // do setup
    if (this->m_KrylovMethod == NULL) {
        ierr = this->SetupKrylovMethod(); CHKERRQ(ierr);
    }

    // invert preconditioner
    ierr = this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);
    ierr = KSPSolve(this->m_KrylovMethod, x, Px); CHKERRQ(ierr);
    ierr = this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief setup 2 level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Setup2LevelPrecond"
PetscErrorCode PrecondReg::Setup2LevelPrecond() {
    PetscErrorCode ierr;
    IntType nl_f,ng_f,nl_c,ng_c,nt,nx_c[3],nx_f[3];
    ScalarType scale,value;
    std::stringstream ss;
    Vec m=NULL,lambda=NULL;
    ScalarType *p_mj=NULL,*p_m=NULL,*p_mjcoarse=NULL,*p_mcoarse=NULL,
                *p_lj=NULL,*p_l=NULL,*p_ljcoarse=NULL,*p_lcoarse=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);


    // check if optimization problem is set up
    ierr = Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_PreProc!=NULL,"null pointer"); CHKERRQ(ierr);

    nt   = this->m_Opt->GetDomainPara().nt;
    nl_f = this->m_Opt->GetDomainPara().nlocal;
    ng_f = this->m_Opt->GetDomainPara().nglobal;

    scale = this->m_Opt->GetKrylovSolverPara().pcgridscale;

    for (int i=0; i < 3; ++i) {
        nx_f[i] = this->m_Opt->GetDomainPara().nx[i];
        value = static_cast<ScalarType>(nx_f[i])/scale;
        nx_c[i] = static_cast<IntType>( std::ceil(value) );
    }

    if (this->m_Opt->GetVerbosity() > 1) {
        ss  << "initializing two-level preconditioner; "
            << "fine level: ("<<nx_f[0]<< ","<<nx_f[1]<< ","<<nx_f[2] << "); "
            << "coarse level: ("<<nx_c[0]<< ","<<nx_c[1]<< ","<<nx_c[2]<< ")";
        ierr = DbgMsg(ss.str()); CHKERRQ(ierr);
        ss.str(std::string()); ss.clear();
    }

    if (this->m_OptProbCoarse==NULL) {

        if (this->m_OptCoarse!=NULL) {
            delete this->m_OptCoarse;
            this->m_OptCoarse=NULL;
        }

        try{ this->m_OptCoarse = new RegOpt(*this->m_Opt); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        for (int i=0; i < 3; ++i) {
            this->m_OptCoarse->SetNumGridPoints(i,nx_c[i]);
        }
        ierr = this->m_OptCoarse->DoSetup(false); CHKERRQ(ierr);

        nl_c = this->m_OptCoarse->GetDomainPara().nlocal;
        ng_c = this->m_OptCoarse->GetDomainPara().nglobal;


        // allocate class for registration
        if (this->m_Opt->GetRegModel() == COMPRESSIBLE) {
            try{ this->m_OptProbCoarse = new OptimalControlRegistration(this->m_OptCoarse); }
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == STOKES) {
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationIC(this->m_OptCoarse); }
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES) {
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationRelaxedIC(this->m_OptCoarse); }
            catch (std::bad_alloc&) {
                ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else{ ierr = ThrowError("registration model not defined"); CHKERRQ(ierr); }


        // create vector fields
        ierr = VecCreate(this->m_WorkScaField1,nl_f,ng_f); CHKERRQ(ierr);
        ierr = VecCreate(this->m_WorkScaField2,nl_f,ng_f); CHKERRQ(ierr);

        try{ this->m_ControlVariable = new VecField(this->m_Opt); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        try{ this->m_IncControlVariable = new VecField(this->m_Opt); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr = VecCreate(this->m_StateVariableCoarse,(nt+1)*nl_c,(nt+1)*ng_c); CHKERRQ(ierr);
        ierr = VecCreate(this->m_AdjointVariableCoarse,(nt+1)*nl_c,(nt+1)*ng_c); CHKERRQ(ierr);

        ierr = VecCreate(this->m_WorkScaFieldCoarse1,nl_c,ng_c); CHKERRQ(ierr);
        ierr = VecCreate(this->m_WorkScaFieldCoarse2,nl_c,ng_c); CHKERRQ(ierr);

        try{ this->m_ControlVariableCoarse = new VecField(this->m_OptCoarse); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        try{ this->m_IncControlVariableCoarse = new VecField(this->m_OptCoarse); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr = VecCreate(this->m_xCoarse,3*nl_c,3*ng_c); CHKERRQ(ierr);
        ierr = VecCreate(this->m_HxCoarse,3*nl_c,3*ng_c); CHKERRQ(ierr);

    }


    // if parameter continuatoin is enabled, parse regularization weight
    if (this->m_Opt->GetParaCont().enabled) {
        this->m_OptCoarse->SetRegularizationWeight(0,this->m_Opt->GetRegNorm().beta[0]);
        this->m_OptCoarse->SetRegularizationWeight(1,this->m_Opt->GetRegNorm().beta[1]);
    }

    // get variables from optimization problem on fine level
    ierr = this->m_OptProb->GetControlVariable(this->m_ControlVariable); CHKERRQ(ierr);
    ierr = this->m_OptProb->GetStateVariable(m); CHKERRQ(ierr);
    ierr = this->m_OptProb->GetAdjointVariable(lambda); CHKERRQ(ierr);

    // restrict control variable
    ierr = this->m_PreProc->Restrict(this->m_ControlVariableCoarse,this->m_ControlVariable,nx_c,nx_f); CHKERRQ(ierr);

    ierr = VecGetArray(m,&p_m); CHKERRQ(ierr);
    ierr = VecGetArray(lambda,&p_l); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_StateVariableCoarse,&p_mcoarse); CHKERRQ(ierr);
    ierr = VecGetArray(this->m_AdjointVariableCoarse,&p_lcoarse); CHKERRQ(ierr);

    nl_c = this->m_OptCoarse->GetDomainPara().nlocal;
    ng_c = this->m_OptCoarse->GetDomainPara().nglobal;

    // apply restriction operator to time series of images
    for (IntType j = 0; j <= nt; ++j) {

        // get time point of state variable on fine grid
        ierr = VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
        try{ std::copy(p_m+j*nl_f, p_m+(j+1)*nl_f, p_mj); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);

        // apply restriction operator to m_j
        ierr = this->m_PreProc->Restrict(&this->m_WorkScaFieldCoarse1,this->m_WorkScaField1,nx_c,nx_f); CHKERRQ(ierr);

        // store restricted state variable
        ierr = VecGetArray(this->m_WorkScaFieldCoarse1, &p_mjcoarse); CHKERRQ(ierr);
        try{ std::copy(p_mjcoarse, p_mjcoarse+nl_c, p_mcoarse+j*nl_c); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_WorkScaFieldCoarse1, &p_mjcoarse); CHKERRQ(ierr);

        // get time point of adjoint variable on fine grid
        ierr = VecGetArray(this->m_WorkScaField2,&p_lj); CHKERRQ(ierr);
        try{ std::copy(p_l+j*nl_f,p_l+(j+1)*nl_f,p_lj); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_WorkScaField2,&p_lj); CHKERRQ(ierr);

        // apply restriction operator
        ierr = this->m_PreProc->Restrict(&this->m_WorkScaFieldCoarse2,this->m_WorkScaField2,nx_c,nx_f); CHKERRQ(ierr);

        // store restricted adjoint variable
        ierr = VecGetArray(this->m_WorkScaFieldCoarse2,&p_ljcoarse); CHKERRQ(ierr);
        try{ std::copy(p_ljcoarse,p_ljcoarse+nl_c,p_lcoarse+j*nl_c); }
        catch(std::exception&) {
            ierr = ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(this->m_WorkScaFieldCoarse2,&p_ljcoarse); CHKERRQ(ierr);

    }

    ierr = VecRestoreArray(this->m_AdjointVariableCoarse,&p_lcoarse); CHKERRQ(ierr);
    ierr = VecRestoreArray(this->m_StateVariableCoarse,&p_mcoarse); CHKERRQ(ierr);
    ierr = VecRestoreArray(lambda,&p_l); CHKERRQ(ierr);
    ierr = VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    // parse variables to optimization problem on coarse level
    ierr = this->m_OptProbCoarse->SetControlVariable(this->m_ControlVariableCoarse); CHKERRQ(ierr);
    ierr = this->m_OptProbCoarse->SetStateVariable(this->m_StateVariableCoarse); CHKERRQ(ierr);
    ierr = this->m_OptProbCoarse->SetAdjointVariable(this->m_AdjointVariableCoarse); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupKrylovMethod"
PetscErrorCode PrecondReg::SetupKrylovMethod() {
    PetscErrorCode ierr = 0;
    PC pc = NULL;
    IntType nl, ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get sizes
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    ierr = Assert(this->m_KrylovMethod == NULL, "expecting null pointer"); CHKERRQ(ierr);
    ierr = KSPCreate(PETSC_COMM_WORLD, &this->m_KrylovMethod); CHKERRQ(ierr);

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver) {
        case CHEB:
        {
            // chebyshev iteration
            ierr = KSPSetType(this->m_KrylovMethod,KSPCHEBYSHEV); CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 7)
            ierr = KSPChebyshevEstEigSetUseRandom(this->m_KrylovMethod,PETSC_TRUE); CHKERRQ(ierr);
#endif
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            ierr = KSPSetType(this->m_KrylovMethod,KSPCG); CHKERRQ(ierr);
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            ierr = KSPSetType(this->m_KrylovMethod,KSPFCG); CHKERRQ(ierr);
            break;
        }
        case GMRES:
        {
            // generalized minimal residual method
            ierr = KSPSetType(this->m_KrylovMethod,KSPGMRES); CHKERRQ(ierr);
            break;
        }
        case FGMRES:
        {
            // flexible generalized minimal residual method
            ierr = KSPSetType(this->m_KrylovMethod,KSPFGMRES); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr = ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    //ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethod,PETSC_FALSE); CHKERRQ(ierr);

    //ierr = KSPSetPostSolve(this->m_KrylovMethod,PostKrylovSolve,this);
    ierr = KSPSetPreSolve(this->m_KrylovMethod,InvertPrecondPreKrylovSolve,this);

    // set up matvec for preconditioner
    if (this->m_MatVec != NULL) {
        ierr = MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr = MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_MatVec); CHKERRQ(ierr);
    ierr = MatShellSetOperation(this->m_MatVec,MATOP_MULT,(void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr = KSPSetOperators(this->m_KrylovMethod,this->m_MatVec,this->m_MatVec);CHKERRQ(ierr);
    ierr = MatSetOption(this->m_MatVec,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
    //ierr = MatSetOption(this->m_MatVec,MAT_SYMMETRIC,PETSC_FALSE); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 1) {
        ierr = KSPMonitorSet(this->m_KrylovMethod,InvertPrecondKrylovMonitor,this,NULL); CHKERRQ(ierr);
    }

    // remove preconditioner
    ierr = KSPGetPC(this->m_KrylovMethod,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr = KSPSetFromOptions(this->m_KrylovMethod); CHKERRQ(ierr);
    ierr = KSPSetUp(this->m_KrylovMethod); CHKERRQ(ierr);


    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for two level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode PrecondReg::HessianMatVec(Vec Hx, Vec x) {
    PetscErrorCode ierr = 0;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // apply hessian (hessian matvec)
    if (this->m_Opt->GetKrylovSolverPara().pctype == TWOLEVEL) {
        ierr = this->HessianMatVecRestrict(Hx, x); CHKERRQ(ierr);
    } else {
        ierr = this->m_OptProb->HessianMatVec(Hx, x); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for two level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVecRestrict"
PetscErrorCode PrecondReg::HessianMatVecRestrict(Vec Hx, Vec x) {
    PetscErrorCode ierr = 0;
    IntType nx_c[3], nx_f[3];
    ScalarType pct, value;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // check if all the necessary pointers have been initialized
    ierr = Assert(this->m_PreProc != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_xCoarse != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_HxCoarse != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_OptProbCoarse != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncControlVariable != NULL, "null pointer"); CHKERRQ(ierr);
    ierr = Assert(this->m_IncControlVariableCoarse != NULL, "null pointer"); CHKERRQ(ierr);

    // allocate vector field
    if (this->m_WorkVecField == NULL) {
        try{ this->m_WorkVecField = new VecField(this->m_Opt); }
        catch (std::bad_alloc&) {
            ierr = reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    pct = 0; // set to zero, cause we search for a max
    for (int i = 0; i < 3; ++i) {
        nx_f[i] = this->m_Opt->GetDomainPara().nx[i];
        nx_c[i] = this->m_OptCoarse->GetDomainPara().nx[i];
        value = static_cast<ScalarType>(nx_c[i]);
        value /= static_cast<ScalarType>(nx_f[i]);
        pct = value > pct ? value : pct;
    }

    // check onput
    ierr = Assert(nx_f[0] >= nx_c[0], "nx_f[0] < nx_c[0]"); CHKERRQ(ierr);
    ierr = Assert(nx_f[1] >= nx_c[1], "nx_f[1] < nx_c[1]"); CHKERRQ(ierr);
    ierr = Assert(nx_f[2] >= nx_c[2], "nx_f[2] < nx_c[2]"); CHKERRQ(ierr);

    // set components
    ierr = this->m_WorkVecField->SetComponents(x); CHKERRQ(ierr);

    // apply low pass filter before we restrict
    // incremental control variable to coarse grid
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_IncControlVariable,
                                                this->m_WorkVecField, pct); CHKERRQ(ierr);

    // apply restriction operator to incremental control variable
    ierr = this->m_PreProc->Restrict(this->m_IncControlVariableCoarse,
                                     this->m_IncControlVariable,
                                     nx_c, nx_f); CHKERRQ(ierr);

    // get the components to interface hessian mat vec
    ierr = this->m_IncControlVariableCoarse->GetComponents(this->m_xCoarse); CHKERRQ(ierr);

    // apply hessian (hessian matvec)
    ierr = this->m_OptProbCoarse->HessianMatVec(this->m_HxCoarse, this->m_xCoarse, false); CHKERRQ(ierr);

    // get components (for interface of hessian matvec)
    ierr = this->m_IncControlVariableCoarse->SetComponents(this->m_HxCoarse); CHKERRQ(ierr);

    // apply prolongation operator
    ierr = this->m_PreProc->Prolong(this->m_IncControlVariable,
                                    this->m_IncControlVariableCoarse,
                                    nx_f, nx_c ); CHKERRQ(ierr);

    // apply low pass filter to output of hessian matvec
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_IncControlVariable,
                                                this->m_IncControlVariable, pct); CHKERRQ(ierr);

    // apply high-pass filter to input
    ierr = this->m_PreProc->ApplyRectFreqFilter(this->m_WorkVecField,
                                                this->m_WorkVecField, pct, false); CHKERRQ(ierr);

    // add up high and low frequency components
    this->m_IncControlVariable->AXPY(1.0, this->m_WorkVecField); CHKERRQ(ierr);

    // parse to output
    ierr = this->m_IncControlVariable->GetComponents(Hx); CHKERRQ(ierr);

    //ierr = this->m_OptProbCoarse->HessianMatVec(Hx,x,false); CHKERRQ(ierr);
    //ierr = this->m_OptProb->HessianMatVec(Hx,x,false); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief this is an interface to compute the eigenvalues needed
 * when considering a chebyshev method to invert the preconditioner;
 * the eigenvalues are estimated using the Lanczo (KSPCG) or
 * Arnoldi (KSPGMRES) process using a random right hand side vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "EstimateEigenValues"
PetscErrorCode PrecondReg::EstimateEigenValues()
{
    PetscErrorCode ierr = 0;
    IntType n, neig, nl, ng;
    KSPConvergedReason reason;
    std::stringstream ss;
    Vec b = NULL, x = NULL;
    ScalarType *re = NULL, *im = NULL, eigmin, eigmax, emin, emax;

    PetscFunctionBegin;

/*
    itermax = 20;

    // get iteration number; if we need to many iterations, we
    // might want to reestimate the eigenvalues
    n = this->m_Opt->GetKrylovSolverPara().iter;
    if ( n > itermax ) {

        ss << "iter="<< n << " > itermax=" << itermax <<" re-estimating eigenvalues";
        ierr = WrngMsg(ss.str()); CHKERRQ(ierr);
        this->m_EigenValuesEstimated=false;

        // reset krylov iterations
        this->m_Opt->SetKrylovIterations(0);

    }
*/
    // if we detect a divergence, we'll have to restimate
    // the eigenvalues
    ierr = KSPGetConvergedReason(this->m_KrylovMethod, &reason); CHKERRQ(ierr);
    if (reason !=  KSP_DIVERGED_ITS) {
        this->m_EigenValuesEstimated = false;
    }

    if (!this->m_EigenValuesEstimated) {

        if (this->m_Opt->GetVerbosity() > 1) {
            ierr = DbgMsg("estimating eigenvalues of hessian"); CHKERRQ(ierr);
        }

        if (this->m_Opt->GetKrylovSolverPara().usepetsceigest) {
            // default interface for chebyshev method; this interface
            // we have to compute the eigenvalues at every iteration
            // using this interface (at least it seems to be like that)
            ierr = KSPChebyshevEstEigSet(this->m_KrylovMethod,PETSC_DECIDE,
                                                              PETSC_DECIDE,
                                                              PETSC_DECIDE,
                                                              PETSC_DECIDE); CHKERRQ(ierr);
        } else {

            // get sizes
            nl = this->m_Opt->GetDomainPara().nlocal;
            ng = this->m_Opt->GetDomainPara().nglobal;

            ierr = VecCreate(x, 3*nl, 3*ng); CHKERRQ(ierr);
            ierr = VecCreate(b, 3*nl, 3*ng); CHKERRQ(ierr);

            // use random right hand side
            if (this->m_RandomNumGen == NULL) {
                ierr = PetscRandomCreate(PetscObjectComm((PetscObject)b),&this->m_RandomNumGen); CHKERRQ(ierr);
            }
            ierr = VecSetRandom(b, this->m_RandomNumGen); CHKERRQ(ierr);

            // do setup
            if (this->m_KrylovMethodEigEst == NULL) {
                ierr = this->SetupKrylovMethodEigEst(); CHKERRQ(ierr);
            }
            ierr = Assert(this->m_KrylovMethodEigEst != NULL, "null pointer"); CHKERRQ(ierr);

            ierr = KSPSolve(this->m_KrylovMethodEigEst, b, x); CHKERRQ(ierr);
            ierr = KSPGetIterationNumber(this->m_KrylovMethodEigEst, &n); CHKERRQ(ierr);

            std::cout << n <<std::endl;
            ierr = PetscMalloc2(n, &re, n, &im); CHKERRQ(ierr);
            ierr = KSPComputeEigenvalues(this->m_KrylovMethodEigEst, n, re, im, &neig); CHKERRQ(ierr);

            eigmin = PETSC_MAX_REAL;
            eigmax = PETSC_MIN_REAL;

            for (IntType i=0; i < neig; ++i) {
                eigmin = PetscMin(eigmin,re[i]);
                eigmax = PetscMax(eigmax,re[i]);
            }

            // clear memory
            ierr = PetscFree2(re, im); CHKERRQ(ierr);
            ierr = this->m_OptProbCoarse->EstimateExtremalHessEigVals(emin, emax); CHKERRQ(ierr);

            ierr = KSPChebyshevSetEigenvalues(this->m_KrylovMethod, eigmax, eigmin); CHKERRQ(ierr);

        }// switch between eigenvalue estimators
    }
    this->m_EigenValuesEstimated = true;

    if (x != NULL) {ierr = VecDestroy(&x); CHKERRQ(ierr);}
    if (b != NULL) {ierr = VecDestroy(&b); CHKERRQ(ierr);}

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief do setup for krylov method to estimate eigenvalues
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupKrylovMethodEigEst"
PetscErrorCode PrecondReg::SetupKrylovMethodEigEst()
{
    PetscErrorCode ierr;
    PC pc=NULL;
    IntType nl,ng;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // get sizes
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;

    // create krylov method
    if (this->m_KrylovMethodEigEst != NULL) {
        ierr = KSPDestroy(&this->m_KrylovMethodEigEst); CHKERRQ(ierr);
        this->m_KrylovMethodEigEst = NULL;
    }
    ierr = KSPCreate(PETSC_COMM_WORLD,&this->m_KrylovMethodEigEst); CHKERRQ(ierr);

    // preconditioned conjugate gradient
    ierr = KSPSetType(this->m_KrylovMethodEigEst,KSPCG); CHKERRQ(ierr);

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr = KSPSetNormType(this->m_KrylovMethodEigEst,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    //ierr = KSPSetNormType(this->m_KrylovMethod,KSP_NORM_PRECONDITIONED); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethodEigEst,PETSC_FALSE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_MatVecEigEst != NULL) {
        ierr = MatDestroy(&this->m_MatVecEigEst); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr = MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_MatVecEigEst); CHKERRQ(ierr);
    ierr = MatShellSetOperation(this->m_MatVecEigEst,MATOP_MULT,(void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr = KSPSetOperators(this->m_KrylovMethodEigEst,this->m_MatVecEigEst,this->m_MatVecEigEst);CHKERRQ(ierr);
    ierr = MatSetOption(this->m_MatVecEigEst,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
    //ierr = MatSetOption(this->m_MatVec,MAT_SYMMETRIC,PETSC_FALSE); CHKERRQ(ierr);

    // remove preconditioner
    ierr = KSPGetPC(this->m_KrylovMethodEigEst,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    ierr = KSPSetTolerances(this->m_KrylovMethodEigEst,1E-12,1E-12,1E+6,10); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(this->m_KrylovMethodEigEst,PETSC_FALSE); CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(this->m_KrylovMethodEigEst,"esteig_"); CHKERRQ(ierr);

    // we are going to estimate eigenvalues with this
    ierr = KSPSetComputeEigenvalues(this->m_KrylovMethodEigEst,PETSC_TRUE); CHKERRQ(ierr);

    // finish
    ierr = KSPSetUp(this->m_KrylovMethodEigEst); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




} // end of namespace




#endif // _PRECONDREG_CPP_
