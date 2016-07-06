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

    this->m_Opt = NULL;
    this->m_OptCoarse = NULL;

    this->m_MatVec = NULL;
    this->m_KrylovMethod = NULL;

    this->m_xCoarse = NULL;
    this->m_HxCoarse = NULL;

    this->m_PreProc = NULL;
    this->m_OptProbCoarse = NULL;

    this->m_ControlVariable = NULL;
    this->m_IncControlVariable = NULL;

    this->m_StateVariableCoarse = NULL;
    this->m_AdjointVariableCoarse = NULL;
    this->m_ControlVariableCoarse = NULL;
    this->m_IncControlVariableCoarse = NULL;

    this->m_WorkScaField1 = NULL; ///< work scalar field
    this->m_WorkScaField2 = NULL; ///< work scalar field
    this->m_WorkScaFieldCoarse1 = NULL; ///< work scalar field (coarse level)
    this->m_WorkScaFieldCoarse2 = NULL; ///< work scalar field (coarse level)

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

    if (this->m_KrylovMethod != NULL){
        ierr=KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }
    if (this->m_MatVec != NULL){
        ierr=MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    if (this->m_xCoarse != NULL){
        ierr=VecDestroy(&this->m_xCoarse); CHKERRQ(ierr);
        this->m_xCoarse = NULL;
    }
    if (this->m_HxCoarse != NULL){
        ierr=VecDestroy(&this->m_HxCoarse); CHKERRQ(ierr);
        this->m_HxCoarse = NULL;
    }
    if (this->m_OptProbCoarse != NULL){
        delete this->m_OptProbCoarse;
        this->m_OptProbCoarse=NULL;
    }
    if (this->m_OptCoarse != NULL){
        delete this->m_OptCoarse;
        this->m_OptCoarse = NULL;
    }
    if (this->m_ControlVariableCoarse != NULL){
        delete this->m_ControlVariableCoarse;
        this->m_ControlVariableCoarse=NULL;
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
PetscErrorCode PrecondReg::SetProblem(PrecondReg::OptProbType* optprob)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(optprob!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_OptProb = optprob;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set the optimization problem (this is a general purpose
 * implementation; the user can set different optimization problems
 * and we can solve them)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetPreProc"
PetscErrorCode PrecondReg::SetPreProc(PreProcReg* preproc)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(preproc!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_PreProc = preproc;

    PetscFunctionReturn(0);
}





/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "MatVec"
PetscErrorCode PrecondReg::MatVec(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // switch case for choice of preconditioner
    switch(this->m_Opt->GetKrylovSolverPara().pctype){
        case NOPC:
        {
            ierr=WrngMsg("no preconditioner used"); CHKERRQ(ierr);
            ierr=VecCopy(x,Px); CHKERRQ(ierr);
            break;
        }
        case INVREG:
        {
            ierr=this->ApplyInvRegPC(Px,x); CHKERRQ(ierr);
            break;
        }
        case TWOLEVEL:
        {
            ierr=this->Apply2LevelPC(Px,x); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner not defined"); CHKERRQ(ierr);
            break;
        }
    }

    // increment counter
    this->m_Opt->IncrementCounter(PCMATVEC);



    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief apply inverse of regularization operator as preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplyInvRegPC"
PetscErrorCode PrecondReg::ApplyInvRegPC(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    // apply inverse regularization operator
    ierr=this->m_OptProb->ApplyInvRegOp(Px,x); CHKERRQ(ierr);

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief applies the preconditioner for the hessian to a vector
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Apply2LevelPC"
PetscErrorCode PrecondReg::Apply2LevelPC(Vec Px, Vec x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);

    // setup preconditioner
    ierr=this->SetUp2LevelPC(); CHKERRQ(ierr);

    // do setup
    if (this->m_KrylovMethod == NULL){
        ierr=this->SetupKrylovMethod(); CHKERRQ(ierr);
    }
    ierr=this->SetTolerancesKrylovMethod(); CHKERRQ(ierr);

    // start timer
    ierr=this->m_Opt->StartTimer(PMVEXEC); CHKERRQ(ierr);

    ierr=KSPSolve(this->m_KrylovMethod,x,Px); CHKERRQ(ierr);

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVEXEC); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief setup 2 level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetUp2LevelPC"
PetscErrorCode PrecondReg::SetUp2LevelPC()
{
    PetscErrorCode ierr;
    IntType nl_f,ng_f,nl_c,ng_c,nt,nx_c[3],nx_f[3];
    ScalarType scale;
    Vec m=NULL,lambda=NULL;
    ScalarType *p_mj=NULL,*p_m=NULL,*p_mjcoarse=NULL,*p_mcoarse=NULL,
                *p_lj=NULL,*p_l=NULL,*p_ljcoarse=NULL,*p_lcoarse=NULL;

    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_OptProb!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_PreProc!=NULL,"null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("setting up two-level preconditioner"); CHKERRQ(ierr);
    }

    // start timer
    ierr=this->m_Opt->StartTimer(PMVSETUP); CHKERRQ(ierr);

    nt = this->m_Opt->GetDomainPara().nt;
    nl_f = this->m_Opt->GetDomainPara().nlocal;
    ng_f = this->m_Opt->GetDomainPara().nglobal;

    scale = 2.0;

    nx_f[0] = this->m_Opt->GetDomainPara().nx[0];
    nx_f[1] = this->m_Opt->GetDomainPara().nx[1];
    nx_f[2] = this->m_Opt->GetDomainPara().nx[2];

    nx_c[0] = static_cast<IntType>( std::ceil( static_cast<ScalarType>(nx_f[0])/scale ) );
    nx_c[1] = static_cast<IntType>( std::ceil( static_cast<ScalarType>(nx_f[1])/scale ) );
    nx_c[2] = static_cast<IntType>( std::ceil( static_cast<ScalarType>(nx_f[2])/scale ) );

    if (this->m_OptProbCoarse==NULL){

        if (this->m_OptCoarse!=NULL){
            delete this->m_OptCoarse;
            this->m_OptCoarse=NULL;
        }

        try{ this->m_OptCoarse = new RegOpt(*this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        this->m_OptCoarse->SetNumGridPoints(0,nx_c[0]);
        this->m_OptCoarse->SetNumGridPoints(1,nx_c[1]);
        this->m_OptCoarse->SetNumGridPoints(2,nx_c[2]);

        ierr=this->m_OptCoarse->DoSetup(false); CHKERRQ(ierr);

        nl_c = this->m_OptCoarse->GetDomainPara().nlocal;
        ng_c = this->m_OptCoarse->GetDomainPara().nglobal;

        // allocate class for registration
        if (this->m_Opt->GetRegModel() == COMPRESSIBLE){
            try{ this->m_OptProbCoarse = new OptimalControlRegistration(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == STOKES){
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationIC(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES){
            try{ this->m_OptProbCoarse = new OptimalControlRegistrationIC(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }
        else{
            try{ this->m_OptProbCoarse = new OptimalControlRegistration(this->m_OptCoarse); }
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr=VecCreate(this->m_WorkScaField1,nl_f,ng_f); CHKERRQ(ierr);
        ierr=VecCreate(this->m_WorkScaField2,nl_f,ng_f); CHKERRQ(ierr);

        try{ this->m_ControlVariable = new VecField(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        try{ this->m_IncControlVariable = new VecField(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr=VecCreate(this->m_StateVariableCoarse,(nt+1)*nl_c,(nt+1)*ng_c); CHKERRQ(ierr);
        ierr=VecCreate(this->m_AdjointVariableCoarse,(nt+1)*nl_c,(nt+1)*ng_c); CHKERRQ(ierr);

        ierr=VecCreate(this->m_WorkScaFieldCoarse1,nl_c,ng_c); CHKERRQ(ierr);
        ierr=VecCreate(this->m_WorkScaFieldCoarse2,nl_c,ng_c); CHKERRQ(ierr);

        try{ this->m_ControlVariableCoarse = new VecField(this->m_OptCoarse); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        try{ this->m_IncControlVariableCoarse = new VecField(this->m_OptCoarse); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }

        ierr=VecCreate(this->m_xCoarse,3*nl_c,3*ng_c); CHKERRQ(ierr);
        ierr=VecCreate(this->m_HxCoarse,3*nl_c,3*ng_c); CHKERRQ(ierr);

    }

    // get variables from optimization problem on fine level
    ierr=this->m_OptProb->GetControlVariable(this->m_ControlVariable); CHKERRQ(ierr);
    ierr=this->m_OptProb->GetStateVariable(m); CHKERRQ(ierr);
    ierr=this->m_OptProb->GetAdjointVariable(lambda); CHKERRQ(ierr);

    ierr=VecGetArray(m,&p_m); CHKERRQ(ierr);
    ierr=VecGetArray(lambda,&p_l); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_StateVariableCoarse,&p_mcoarse); CHKERRQ(ierr);
    ierr=VecGetArray(this->m_AdjointVariableCoarse,&p_lcoarse); CHKERRQ(ierr);

    // apply restriction operator
    for (IntType j = 0; j <= nt; ++j){

        // get current time point
        ierr=VecGetArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);
        try{ std::copy(p_m+j*nl_f,p_m+(j+1)*nl_f,p_mj); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_WorkScaField1,&p_mj); CHKERRQ(ierr);

        // get current time point
        ierr=VecGetArray(this->m_WorkScaField2,&p_lj); CHKERRQ(ierr);
        try{ std::copy(p_l+j*nl_f,p_l+(j+1)*nl_f,p_lj); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_WorkScaField2,&p_lj); CHKERRQ(ierr);

        ierr=this->m_PreProc->Restrict(&this->m_WorkScaFieldCoarse1,this->m_WorkScaField1,nx_c,nx_f); CHKERRQ(ierr);
        ierr=this->m_PreProc->Restrict(&this->m_WorkScaFieldCoarse2,this->m_WorkScaField2,nx_c,nx_f); CHKERRQ(ierr);

        // get current time point
        ierr=VecGetArray(this->m_WorkScaFieldCoarse1,&p_mjcoarse); CHKERRQ(ierr);
        try{ std::copy(p_mjcoarse,p_mjcoarse+nl_c,p_mcoarse+j*nl_c); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_WorkScaFieldCoarse1,&p_mjcoarse); CHKERRQ(ierr);

        // get current time point
        ierr=VecGetArray(this->m_WorkScaFieldCoarse2,&p_ljcoarse); CHKERRQ(ierr);
        try{ std::copy(p_ljcoarse,p_ljcoarse+nl_c,p_lcoarse+j*nl_c); }
        catch(std::exception&){
            ierr=ThrowError("copy failed"); CHKERRQ(ierr);
        }
        ierr=VecRestoreArray(this->m_WorkScaFieldCoarse2,&p_ljcoarse); CHKERRQ(ierr);
    }

    ierr=VecRestoreArray(this->m_AdjointVariableCoarse,&p_lcoarse); CHKERRQ(ierr);
    ierr=VecRestoreArray(this->m_StateVariableCoarse,&p_mcoarse); CHKERRQ(ierr);
    ierr=VecRestoreArray(lambda,&p_l); CHKERRQ(ierr);
    ierr=VecRestoreArray(m,&p_m); CHKERRQ(ierr);

    // restrict control variable
    ierr=this->m_PreProc->Restrict(this->m_ControlVariableCoarse,this->m_ControlVariable,nx_c,nx_f); CHKERRQ(ierr);

//    ierr=this->m_ControlVariableCoarse->Copy(this->m_ControlVariable); CHKERRQ(ierr);
//    ierr=VecCopy(m,this->m_StateVariableCoarse); CHKERRQ(ierr);
//    ierr=VecCopy(lambda,this->m_AdjointVariableCoarse); CHKERRQ(ierr);

    // parse variables to optimization problem on coarse level
    ierr=this->m_OptProbCoarse->SetControlVariable(this->m_ControlVariableCoarse); CHKERRQ(ierr);
    ierr=this->m_OptProbCoarse->SetStateVariable(this->m_StateVariableCoarse); CHKERRQ(ierr);
    ierr=this->m_OptProbCoarse->SetAdjointVariable(this->m_AdjointVariableCoarse); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("setup of two-level preconditioner done"); CHKERRQ(ierr);
    }

    // stop timer
    ierr=this->m_Opt->StopTimer(PMVSETUP); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupKrylovMethod"
PetscErrorCode PrecondReg::SetupKrylovMethod()
{
    PetscErrorCode ierr;
    PC pc=NULL;
    IntType nl,ng;

    PetscFunctionBegin;

    // get sizes
    nl = this->m_Opt->GetDomainPara().nlocal;
    ng = this->m_Opt->GetDomainPara().nglobal;


    // create krylov method
    if (this->m_KrylovMethod != NULL){
        ierr=KSPDestroy(&this->m_KrylovMethod); CHKERRQ(ierr);
        this->m_KrylovMethod = NULL;
    }
    ierr=KSPCreate(PETSC_COMM_WORLD,&this->m_KrylovMethod); CHKERRQ(ierr);


    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            ierr=KSPSetType(this->m_KrylovMethod,KSPCHEBYSHEV); CHKERRQ(ierr);
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPCG); CHKERRQ(ierr);
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            ierr=KSPSetType(this->m_KrylovMethod,KSPFCG); CHKERRQ(ierr);
            break;
        }
        case GMRES:
        {
            // GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPGMRES); CHKERRQ(ierr);
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            ierr=KSPSetType(this->m_KrylovMethod,KSPFGMRES); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }

    //KSP_NORM_UNPRECONDITIONED unpreconditioned norm: ||b-Ax||_2)
    //KSP_NORM_PRECONDITIONED   preconditioned norm: ||P(b-Ax)||_2)
    //KSP_NORM_NATURAL          natural norm: sqrt((b-A*x)*P*(b-A*x))
    ierr=KSPSetNormType(this->m_KrylovMethod,KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
    ierr=KSPSetInitialGuessNonzero(this->m_KrylovMethod,PETSC_TRUE); CHKERRQ(ierr);

    // set up matvec for preconditioner
    if (this->m_MatVec != NULL){
        ierr=MatDestroy(&this->m_MatVec); CHKERRQ(ierr);
        this->m_MatVec = NULL;
    }

    ierr=MatCreateShell(PETSC_COMM_WORLD,3*nl,3*nl,3*ng,3*ng,this,&this->m_MatVec); CHKERRQ(ierr);
    ierr=MatShellSetOperation(this->m_MatVec,MATOP_MULT,(void(*)(void))InvertPrecondMatVec); CHKERRQ(ierr);
    ierr=KSPSetOperators(this->m_KrylovMethod,this->m_MatVec,this->m_MatVec);CHKERRQ(ierr);
    ierr=KSPMonitorSet(this->m_KrylovMethod,InvertPrecondKrylovMonitor,this,PETSC_NULL); CHKERRQ(ierr);

    // remove preconditioner
    ierr=KSPGetPC(this->m_KrylovMethod,&pc); CHKERRQ(ierr);
    ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr); ///< set no preconditioner

    // finish
    ierr=KSPSetUp(this->m_KrylovMethod); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set the tolerances for krylov method
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetTolerancesKrylovMethod"
PetscErrorCode PrecondReg::SetTolerancesKrylovMethod()
{
    PetscErrorCode ierr;
    ScalarType reltol,abstol,divtol,scale;
    IntType maxit;

    PetscFunctionBegin;

    divtol = 1E+06;
    abstol = 1E-16;
    reltol = 1E-16;
    maxit  = 1000;

    // check for null pointer
    ierr=Assert(this->m_KrylovMethod!=NULL,"null pointer"); CHKERRQ(ierr);

    scale = this->m_Opt->GetKrylovSolverPara().pcsolvertol;

    switch (this->m_Opt->GetKrylovSolverPara().pcsolver){
        case CHEB:
        {
            // chebyshev iteration
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        case PCG:
        {
            // preconditioned conjugate gradient
            reltol = scale*this->m_Opt->GetKrylovSolverPara().reltol;
            break;
        }
        case FCG:
        {
            // flexible conjugate gradient
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        case GMRES:
        {
            // GMRES
            reltol = scale*this->m_Opt->GetKrylovSolverPara().reltol;
            break;
        }
        case FGMRES:
        {
            // flexible GMRES
            maxit = this->m_Opt->GetKrylovSolverPara().pcsolvermaxit;
            break;
        }
        default:
        {
            ierr=ThrowError("preconditioner solver not defined"); CHKERRQ(ierr);
            break;
        }
    }
    reltol = std::max(reltol,1E-16); // make sure tolerance is non-zero
    reltol = std::min(reltol,0.5); // make sure tolerance smaller than 0.25

    // set tolerances
    ierr=KSPSetTolerances(this->m_KrylovMethod,reltol,abstol,divtol,maxit); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for two level preconditioner
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "HessianMatVec"
PetscErrorCode PrecondReg::HessianMatVec(Vec Hx, Vec x)
{
    PetscErrorCode ierr;
    IntType nx_c[3],nx_f[3];
    PetscFunctionBegin;

    // check if optimization problem is set up
    ierr=Assert(this->m_xCoarse!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_HxCoarse!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_PreProc!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_OptProbCoarse!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncControlVariable!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_IncControlVariableCoarse!=NULL,"null pointer"); CHKERRQ(ierr);

    // get number of grid points
    nx_c[0] = this->m_OptCoarse->GetDomainPara().nx[0];
    nx_c[1] = this->m_OptCoarse->GetDomainPara().nx[1];
    nx_c[2] = this->m_OptCoarse->GetDomainPara().nx[2];

    nx_f[0] = this->m_Opt->GetDomainPara().nx[0];
    nx_f[1] = this->m_Opt->GetDomainPara().nx[1];
    nx_f[2] = this->m_Opt->GetDomainPara().nx[2];

    ierr=Assert(nx_f[0] >= nx_c[0],"nx_f[0] < nx_c[0]"); CHKERRQ(ierr);
    ierr=Assert(nx_f[1] >= nx_c[1],"nx_f[1] < nx_c[1]"); CHKERRQ(ierr);
    ierr=Assert(nx_f[2] >= nx_c[2],"nx_f[2] < nx_c[2]"); CHKERRQ(ierr);

    // set components
    ierr=this->m_IncControlVariable->SetComponents(x); CHKERRQ(ierr);

    // apply restriction operator
    ierr=this->m_PreProc->Restrict(this->m_IncControlVariableCoarse,this->m_IncControlVariable,nx_c,nx_f); CHKERRQ(ierr);
    ierr=this->m_IncControlVariableCoarse->GetComponents(this->m_xCoarse); CHKERRQ(ierr);

    // apply hessian (hessian matvec)
    ierr=this->m_OptProbCoarse->HessianMatVec(this->m_HxCoarse,this->m_xCoarse); CHKERRQ(ierr);

    // get components (for interface of hessian matvec)
    ierr=this->m_IncControlVariableCoarse->SetComponents(this->m_HxCoarse); CHKERRQ(ierr);

    // apply prolongation operator
    ierr=this->m_PreProc->Prolong(this->m_IncControlVariable,this->m_IncControlVariableCoarse,nx_f,nx_c); CHKERRQ(ierr);
    ierr=this->m_IncControlVariable->GetComponents(Hx); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




} // end of namespace




#endif // _PRECONDREG_CPP_
