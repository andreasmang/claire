/*
 * @file RegistrationInterface.cpp
 *
 * @author Andreas Mang
 */


#ifndef _REGISTRATIONINTERFACE_CPP_
#define _REGISTRATIONINTERFACE_CPP_

#include "RegistrationInterface.hpp"


namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegistrationInterface"
RegistrationInterface::RegistrationInterface()
{
    this->Initialize();
}




/********************************************************************
 * @brief default destructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~RegistrationInterface"
RegistrationInterface::~RegistrationInterface(void)
{
    this->ClearMemory();
}




/********************************************************************
 * @brief constructor
 * @param opt base class for registration options and arguments
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RegistrationInterface"
RegistrationInterface::RegistrationInterface(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief init variables
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode RegistrationInterface::Initialize(void)
{
    PetscFunctionBegin;

    this->m_Opt=NULL;
    this->m_PreProc=NULL;
    this->m_ReadWrite=NULL;
    this->m_Optimizer=NULL;
    this->m_Precond=NULL;
    this->m_RegProblem=NULL;
    this->m_ReferencePyramid=NULL;
    this->m_TemplatePyramid=NULL;

    this->m_Solution=NULL;
    this->m_ReferenceImage=NULL;
    this->m_TemplateImage=NULL;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief clean up
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode RegistrationInterface::ClearMemory(void)
{
    PetscErrorCode ierr=0;
    PetscFunctionBegin;

    // delete class for registration problem
    if (this->m_RegProblem != NULL){
        delete this->m_RegProblem;
        this->m_RegProblem = NULL;
    }

    // delete class for optimizer
    if (this->m_Optimizer != NULL){
        delete this->m_Optimizer;
        this->m_Optimizer = NULL;
    }

    // delete class for pre-processing
    if (this->m_PreProc != NULL){
        delete this->m_PreProc;
        this->m_PreProc = NULL;
    }

    // delete class for preconditioner
    if (this->m_Precond != NULL){
        delete this->m_Precond;
        this->m_Precond = NULL;
    }

    if (this->m_Solution != NULL){
        delete this->m_Solution;
        this->m_Solution=NULL;
    }

    if (this->m_ReferencePyramid != NULL){
        delete this->m_ReferencePyramid;
        this->m_ReferencePyramid=NULL;
    }

    if (this->m_TemplatePyramid != NULL){
        delete this->m_TemplatePyramid;
        this->m_TemplatePyramid=NULL;
    }

    // if we did not read the images, we can
    // destroy the containers here
    if (!this->m_Opt->GetReadWriteFlags().readfiles){

        // delete reference image
        if (this->m_ReferenceImage != NULL){
            ierr=VecDestroy(&this->m_ReferenceImage); CHKERRQ(ierr);
            this->m_ReferenceImage = NULL;
        }

        // delete template image
        if (this->m_TemplateImage != NULL){
            ierr=VecDestroy(&this->m_TemplateImage); CHKERRQ(ierr);
            this->m_TemplateImage = NULL;
        }

    }

    PetscFunctionReturn(ierr);
}




/********************************************************************
 * @brief set initial guess
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode RegistrationInterface::SetInitialGuess(VecField* x)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // the input better is not zero
    ierr=Assert(x!=NULL,"null pointer"); CHKERRQ(ierr);

    // if we have not setup initial guess, do so
    if (this->m_Solution == NULL){
        try{this->m_Solution = new VecField(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    ierr=this->m_Solution->Copy(x); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReadWrite"
PetscErrorCode RegistrationInterface::SetReadWrite(ReadWriteReg* rw)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(rw!=NULL,"null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite = rw;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set reference image (i.e., the fixed image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReferenceImage"
PetscErrorCode RegistrationInterface::SetReferenceImage(Vec mR)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(mR!=NULL,"reference image is null"); CHKERRQ(ierr);
    ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);

    this->m_ReferenceImage = mR;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set template image (i.e., the template image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetTemplateImage"
PetscErrorCode RegistrationInterface::SetTemplateImage(Vec mT)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(mT!=NULL,"template image is null"); CHKERRQ(ierr);
    ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

    this->m_TemplateImage = mT;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set read/write object
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "DispLevelMsg"
PetscErrorCode RegistrationInterface::DispLevelMsg(std::string msg, int rank)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg(msg); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set up the registration problem and optimizer
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupSolver"
PetscErrorCode RegistrationInterface::SetupSolver()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // reset optimization problem
    if (this->m_Optimizer != NULL){
        delete this->m_Optimizer; this->m_Optimizer = NULL;
    }

    // allocate class for io
    try{ this->m_Optimizer = new OptimizerType(this->m_Opt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // set up optimization problem
    ierr=this->SetupRegProblem(); CHKERRQ(ierr);

    if(this->m_PreProc!=NULL){
        delete this->m_PreProc; this->m_PreProc=NULL;
    }
    try{this->m_PreProc = new PreProcReg(this->m_Opt);}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    if(this->m_Opt->GetKrylovSolverPara().pctype != NOPC){

        if(this->m_Precond!=NULL){
            delete this->m_Precond; this->m_Precond=NULL;
        }
        // allocate preconditioner
        try{ this->m_Precond = new PrecondReg(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        ierr=this->m_Precond->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
        ierr=this->m_Precond->SetProblem(this->m_RegProblem); CHKERRQ(ierr);
        ierr=this->m_Optimizer->SetPreconditioner(this->m_Precond); CHKERRQ(ierr);

    }

    // set up initial condition
    if (this->m_Solution==NULL){

        try{ this->m_Solution = new VecField(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
        this->m_Solution->SetValue(0.0); CHKERRQ(ierr);

    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief set up the registration problem, which essentially is
 * equivalent to allocating the class
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupRegProblem"
PetscErrorCode RegistrationInterface::SetupRegProblem()
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // reset registration problem
    if (this->m_RegProblem != NULL){
        delete this->m_RegProblem; this->m_RegProblem = NULL;
    }

    // allocate class for registration
    if (this->m_Opt->GetRegModel() == COMPRESSIBLE){
        try{ this->m_RegProblem = new OptimalControlRegistration(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    else if (this->m_Opt->GetRegModel() == STOKES){
        try{ this->m_RegProblem = new OptimalControlRegistrationIC(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    else if (this->m_Opt->GetRegModel() == RELAXEDSTOKES){
        try{ this->m_RegProblem = new OptimalControlRegistrationRelaxedIC(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    else{ ierr=ThrowError("registration modle not available"); CHKERRQ(ierr); }

    ierr=Assert(this->m_ReadWrite!=NULL,"read/write is null"); CHKERRQ(ierr);
    ierr=this->m_RegProblem->SetReadWrite(this->m_ReadWrite); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief main function to call in order to solve the optimization
 * problem
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Run"
PetscErrorCode RegistrationInterface::Run()
{
    PetscErrorCode ierr;
    IntType nxmax,nx;
    std::stringstream ss;
    int rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr=Msg("starting optimization"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=PetscPrintf(PETSC_COMM_WORLD," %s  %-20s %-20s %-20s %-20s %-20s\n",
                                    "iter","objective (rel)","mismatch (rel)",
                                    "||gradient||_2,rel","||gradient||_2","step"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;


    // switch between solvers we have to solve optimization problem
    if( this->m_Opt->GetParaCont().enabled ){

        // run different flavours of parameter continuation
        ierr=this->RunSolverRegParaCont(); CHKERRQ(ierr);

    }
    else if( this->m_Opt->GetScaleContPara().enabled ){

        // run different flavours of parameter continuation
        ierr=this->RunSolverScaleCont(); CHKERRQ(ierr);

    }
    else if( this->m_Opt->GetGridContPara().enabled ){

        nxmax = PETSC_MIN_INT;
        for (int i=0; i < 3; ++i){
            nx = this->m_Opt->GetDomainPara().nx[i];
            nxmax = nx > nxmax ? nx : nxmax;
        }

        // run grid continuation
        if (nxmax >= 32){ ierr=this->RunSolverGridCont(); CHKERRQ(ierr); }
        else{

            ss << "max(nx) = " << nxmax << " too small for grid continuation; switching to default solver";
            ierr=WrngMsg(ss.str()); CHKERRQ(ierr);
            ss.str(std::string()); ss.clear();

            ierr=this->RunSolver(); CHKERRQ(ierr);
        }
    }
    else{ ierr=this->RunSolver(); CHKERRQ(ierr); }

    ierr=this->DispLevelMsg("optimization done",rank); CHKERRQ(ierr);

    ierr=this->Finalize(); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run single level solver (no grid, scale, or parameter
 * continuation is performed)
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolver"
PetscErrorCode RegistrationInterface::RunSolver()
{
    PetscErrorCode ierr;
    Vec mT=NULL,mR=NULL,x=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // do the setup
    ierr=this->SetupSolver(); CHKERRQ(ierr);

    ierr=Assert(this->m_RegProblem!= NULL, "registration problem is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Optimizer!= NULL, "optimizer is null"); CHKERRQ(ierr);

    // presmoothing, if necessary
    if (this->m_Opt->GetReadWriteFlags().readfiles){

        ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
        ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

        // allocate
        ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);
        ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);

        if (this->m_Opt->GetRegFlags().smoothingenabled){

            if(this->m_PreProc==NULL){
                try{this->m_PreProc = new PreProcReg(this->m_Opt);}
                catch (std::bad_alloc&){
                    ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }

            ierr=this->m_PreProc->ApplySmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
            ierr=this->m_PreProc->ApplySmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

        }
        else{
            ierr=VecCopy(this->m_ReferenceImage,mR); CHKERRQ(ierr);
            ierr=VecCopy(this->m_TemplateImage,mT); CHKERRQ(ierr);
        }

        // rescale images
        ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
        ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

        ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{
        // set up synthetic test problem
        ierr=this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage,this->m_TemplateImage); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);

    }

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // initialize registration problem (evaluate objective and gradient
    // for zero velocity field)
    ierr=this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    // init solver
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // run the optimization
    ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

    // finalize the registration
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vectors
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief main function to run the parameter continuation
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverRegParaCont"
PetscErrorCode RegistrationInterface::RunSolverRegParaCont()
{
    PetscErrorCode ierr;
    Vec mT=NULL,mR=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // do the setup
    ierr=this->SetupSolver(); CHKERRQ(ierr);

    // check if setup was complete
    ierr=Assert(this->m_RegProblem!= NULL, "registration problem is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Optimizer!= NULL, "optimizer is null"); CHKERRQ(ierr);

    // presmoothing, if necessary
    if (this->m_Opt->GetReadWriteFlags().readfiles){

        ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
        ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

        // allocate
        ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);
        ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);

        if (this->m_Opt->GetRegFlags().smoothingenabled){

            if(this->m_PreProc==NULL){
                try{this->m_PreProc = new PreProcReg(this->m_Opt);}
                catch (std::bad_alloc&){
                    ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
                }
            }

            ierr=this->m_PreProc->ApplySmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
            ierr=this->m_PreProc->ApplySmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

        }
        else{
            ierr=VecCopy(this->m_ReferenceImage,mR); CHKERRQ(ierr);
            ierr=VecCopy(this->m_TemplateImage,mT); CHKERRQ(ierr);
        }

        // rescale images
        ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
        ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

        ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{
        // set up synthetic test problem
        ierr=this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage,this->m_TemplateImage); CHKERRQ(ierr);

        ierr=this->m_RegProblem->SetReferenceImage(this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(this->m_TemplateImage); CHKERRQ(ierr);

    }

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // initialize registration problem (evaluate objective and gradient
    // for zero velocity field)
    ierr=this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

    // switch between the different strategies for
    // doing the parameter continuation (default one
    // is binary search)
    switch (this->m_Opt->GetParaCont().strategy){
        case PCONTBINSEARCH:
        {
            ierr=this->RunSolverRegParaContBinarySearch(); CHKERRQ(ierr);
            break;
        }
        case PCONTREDUCESEARCH:
        {
            ierr=this->RunSolverRegParaContReductSearch(); CHKERRQ(ierr);
            break;
        }
        case PCONTINUATION:
        {
            ierr=this->RunSolverRegParaContReduction(); CHKERRQ(ierr);
            break;
        }
        default:
        {
            ierr=ThrowError("parameter continuation strategy not valid"); CHKERRQ(ierr);
            break;
        }
    }

    // destroy vector
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run the optimizer; we search for an optimal
 * regularization weight using a binary search; we reduce/lift the
 * regularization parameter until we found a deformation map that
 * is diffeomorphic and results in a map that is close to the bound
 * on jacobian set by user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverRegParaContBinarySearch"
PetscErrorCode RegistrationInterface::RunSolverRegParaContBinarySearch()
{
    PetscErrorCode ierr;
    int maxsteps,level,rank;
    bool stop,boundreached,converged;
    std::ofstream logwriter;
    std::stringstream ss;
    std::string filename;
    ScalarType beta,betamin,betascale,dbetascale,
                betastar,betahat,dbeta,dbetamin;
    Vec x;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(this->m_Optimizer!=NULL,"optimizer is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!=NULL,"registration problem is null"); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    betamin   = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetParaCont().betascale;
    maxsteps  = this->m_Opt->GetParaCont().maxsteps;

    ierr=Assert(betascale < 1.0 && betascale > 0.0,"scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0 && betamin < 1.0,"lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // initialize parameters
    beta = 1.0;
    betastar = beta;

    // reduce regularization parameter by one order of magnitude until
    // we hit tolerance
    stop=false; level = 0;
    while(level < maxsteps){

        this->m_Opt->SetRegularizationWeight(0,beta);
        this->m_Opt->SetRegularizationWeight(1,beta);
        //this->m_Opt->InitialGradNormSet(false);


        ss << std::scientific << std::setw(3)
            << "level "<< level <<" ( betav="<<beta
            <<"; betav*="<<betastar<<" )";
        ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
        ss.str( std::string() ); ss.clear();

        // set initial guess for current level
        ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);
        ierr=this->m_Optimizer->GetSolutionStatus(converged); CHKERRQ(ierr);
        if (!converged) break;

        // get the solution
        ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

        // check bounds on jacobian
        ierr=this->m_RegProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember regularization parameter
        betastar = beta;

        // if we got here, the solution is valid
        ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

        // reduce beta
        beta *= betascale;

        // if regularization parameter is smaller than
        // lower bound, let's stop this
        if (beta < betamin){

            if (this->m_Opt->GetVerbosity() > 0){
                ss << std::scientific
                   <<"regularization parameter smaller than lower bound (betav="
                   <<beta<<" < " << betamin <<"=betavmin)";
                ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str( std::string() ); ss.clear();
            }
            break;
        }

        ++level;

    } // until we hit the tolerance

    // now start binary search
    stop=false;

    // get scale for delta beta; this parameter determines how
    // accurate we solve (how many digits) with respect to the
    // order of magnitude of magnitude of the regularization
    // parameter)
    dbetascale = this->m_Opt->GetParaCont().dbetascale;
    ierr=Assert(dbetascale<1.0 && dbetascale>0.0,"scale for delta betav not in (0,1)"); CHKERRQ(ierr);

    //update beta
    dbetamin = dbetascale*betastar;
    betahat  = betascale*betastar;
    dbeta    = (betastar-betahat)/2.0;
    beta     = betastar-dbeta;

    ++level;

    while(!stop){

        // set regularization parameter
        this->m_Opt->SetRegularizationWeight(0,beta);
        this->m_Opt->SetRegularizationWeight(1,beta);

        // display regularization parameter to user
        ss<<std::setw(3)<<"level "<<level<<" ( betav="<<beta<<"; betav*="<<betastar<<" )";
        ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
        ss.str( std::string() ); ss.clear();

        // set initial guess for current level
        ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        // if we did not converge, beta is too small, also
        ierr=this->m_Optimizer->GetSolutionStatus(converged); CHKERRQ(ierr);
        if (converged){

            // get the solution
            ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

            // check bounds on jacobian
            boundreached=false;
            ierr=this->m_RegProblem->CheckBounds(x,boundreached); CHKERRQ(ierr);

            // if bound is reached, the lower bound is now beta
            // if not, beta is our new best estimate
            if (boundreached){ betahat = beta; }
            else{

                betastar = beta; // new best estimate

                // if we got here, the solution is valid
                ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

            }
        }
        else{
            ierr=WrngMsg("solver did not converge"); CHKERRQ(ierr);
            betahat = beta;
        };

        // increase or reduce beta
        dbeta = (betastar - betahat)/2.0;
        beta  = betastar - dbeta;
        if (fabs(dbeta) < dbetamin){
            stop = true;
            if (this->m_Opt->GetVerbosity() > 0){
                ss  << std::setw(3)<<"update for beta to small ( "
                    << fabs(dbeta) << " < " << dbetamin << " )";
                ierr=DbgMsg(ss.str());
                ss.str( std::string() ); ss.clear();
            }
        }

        ++level;

    }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss <<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();


    if (rank == 0){
        filename  = this->m_Opt->GetReadWriteFlags().xfolder;
        filename += "parameter-continuation-estimated-beta.log";
        // create output file or append to output file
        logwriter.open(filename.c_str(), std::ofstream::out | std::ofstream::app );
        ierr=Assert(logwriter.is_open(),"could not open file for writing"); CHKERRQ(ierr);
        ss  << std::scientific <<  "betav " << std::setw(3) << std::right << betastar;
        logwriter << ss.str() << std::endl;
        ss.str( std::string() ); ss.clear();
    }

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief we solves the optimization problem by simply reducing
 * the regularization parameter until the mapping becomes
 * non-diffeomorphic/breaches the user defined bound; stored
 * velocity field (solution) is last iterate that resulted in
 * diffeomorphic deformation map (as judged by the determinant
 * of the deformation gradient)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverRegParaContReductSearch"
PetscErrorCode RegistrationInterface::RunSolverRegParaContReductSearch()
{
    PetscErrorCode ierr;
    std::stringstream ss;
    ScalarType beta,betamin,betastar,betascale;
    Vec x;
    int maxsteps,rank,level;
    bool stop;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    betamin   = this->m_Opt->GetBetaMinParaCont();
    maxsteps  = this->m_Opt->GetParaCont().maxsteps;
    betascale = this->m_Opt->GetParaCont().betascale;

    ierr=Assert(betascale < 1.0 && betascale > 0.0,"scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0 && betamin < 1.0,"lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // set initial regularization weight
    beta = 1.0;
    betastar = beta;

    // reduce regularization parameter by one order of magnitude until
    // we hit user defined tolerances (which either is a lower bound
    // on the regularization parameter or a lower bound on the
    // determinant of the deformation gradient)
    level = 0;
    while(level < maxsteps){

        // set regularization weight
        this->m_Opt->SetRegularizationWeight(0,beta);
        this->m_Opt->SetRegularizationWeight(1,beta);

        // display message to user
        ss << std::scientific<<std::setw(3)<<"level "<<level<<" (beta="<<beta<<")";
        ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
        ss.str( std::string() ); ss.clear();

        // set initial guess for current level
        ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        // get the solution
        ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

        // check bounds on jacobian
        stop=false;
        ierr=this->m_RegProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember best estimate
        betastar = beta;

        // if we got here, the solution is valid
        ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

        beta *= betascale; // reduce beta

        // if the regularization parameter is smaller than
        // the lower bound, we're done
        if (beta < betamin){

            if (this->m_Opt->GetVerbosity() > 0){
                ss <<"regularization parameter smaller than lower bound (betav="
                   <<beta<<" < " << betamin <<"=betavmin)";
                ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str( std::string() ); ss.clear();
            }

            break;

        }

        ++level;

    } // parameter reduction

    ss<<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
    ss.str( std::string() ); ss.clear();

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief we solves the optimization problem by simply reducing
 * the regularization parameter until we have reached the
 * target regularization weight set by the user
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverRegParaContReduction"
PetscErrorCode RegistrationInterface::RunSolverRegParaContReduction()
{
    PetscErrorCode ierr;
    std::stringstream ss;
    ScalarType beta,betastar;
    Vec x;
    int level,rank;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get target regularization weight
    betastar=this->m_Opt->GetRegNorm().beta[0];
    ierr=Assert(betastar>0.0 && betastar<1.0,"target beta not in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // reduce regularization parameter
    level = 0; beta=1.0;
    while(beta > betastar){

        // set regularization weight
        this->m_Opt->SetRegularizationWeight(0,beta);
        this->m_Opt->SetRegularizationWeight(1,beta);

        // display message to user
        ss << std::scientific<<std::setw(3)<<"level "<<level<<" (beta="<<beta<<"; beta*="<<betastar<<")";
        ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
        ss.str( std::string() ); ss.clear();

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        beta /= static_cast<ScalarType>(10); // reduce beta
        ++level;

    } // parameter reduction

    beta = betastar;

    // set regularization weight
    this->m_Opt->SetRegularizationWeight(0,beta);
    this->m_Opt->SetRegularizationWeight(1,beta);

    // display message to user
    ss << std::scientific << std::setw(3)
        <<"level "<< level <<" (beta="
        <<beta<<"; beta*="<<betastar<<")";
    ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
    ss.str( std::string() ); ss.clear();

    // solve optimization problem for user defined regularization parameter
    ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run solver using a scale continuation scheme; that is, we
 * will successively reduce the smoothing of the images to be
 * registered to get to finer and finer scales; this is supposed to
 * reduce the non-linearity in the problem
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverScaleCont"
PetscErrorCode RegistrationInterface::RunSolverScaleCont()
{
    PetscErrorCode ierr;
    Vec mT=NULL,mR=NULL,x=NULL;
    std::stringstream ss;
    int level,maxlevel=6,rank;
    bool solve;
    ScalarType sigma[3],nxhalf;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // set up preprocessing
    if(this->m_PreProc==NULL){
        try{this->m_PreProc = new PreProcReg(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // do the setup
    ierr=this->SetupSolver(); CHKERRQ(ierr);

    // check if everything has been set up correctly
    ierr=Assert(this->m_Optimizer!= NULL,"optimizer is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!= NULL,"registration problem is null"); CHKERRQ(ierr);

    // set up synthetic problem if we did not read images
    if (!this->m_Opt->GetReadWriteFlags().readfiles){
        ierr=this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage,this->m_TemplateImage); CHKERRQ(ierr);
    }

    // check if images have been set
    ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

    // allocate local images
    ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);
    ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);

    // set images
    ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    level=0;
    while (level < maxlevel){

        solve=true;
        for (int i=0; i < 3; ++i){

            // get and set sigma for current level
            sigma[i] = this->m_Opt->GetScaleContPara().sigma[i][level];
            this->m_Opt->SetSigma(i,sigma[i]);

            // if sigma is bigger than half of the grid size, don't compute
            nxhalf = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i])/2.0;
            if ( nxhalf <= sigma[i] ) solve=false;

        }

        // solve problem
        if (solve){

            ierr=this->m_PreProc->ApplySmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
            ierr=this->m_PreProc->ApplySmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

            // rescale images
            ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
            ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

            // display message to user
            ss << std::scientific << std::setw(3)
                <<"level "<<level<<" sigma=("
                <<sigma[0]<<","<<sigma[1]<<","<<sigma[2]<<")";
            ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();

            // compute gradient, distance measure, and initial objective
            // value for zero velocity field, but updated images
            ierr=this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

            // run the optimization
            ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        }
        else{

            ss << std::scientific << std::setw(3)
                << "skipping level "<< level <<" sigma=("
                <<sigma[0]<<","<<sigma[1]<<","<<sigma[2]<<")";
            ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();

        }

        ++level; // increment counter

    }

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=this->m_Solution->SetComponents(x); CHKERRQ(ierr);

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vector
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief run solver using a grid continuation scheme; that is, we
 * will successively increase the number of grid points of the
 * template and reference image
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunSolverGridCont"
PetscErrorCode RegistrationInterface::RunSolverGridCont()
{
    PetscErrorCode ierr;
    int rank,level,nlevels;
    std::stringstream ss;
    IntType nx[3],nl,ng;
    Vec mT=NULL,mR=NULL,xstar=NULL;
    VecField *v=NULL;
    ScalarType greltol,tolscale=1E1;
    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // set up preprocessing
    if(this->m_PreProc==NULL){
        try{this->m_PreProc = new PreProcReg(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    this->m_PreProc->ResetGridChangeOps(true);

    if (!this->m_Opt->GetReadWriteFlags().readfiles){

        // do the setup
        ierr=this->SetupSolver(); CHKERRQ(ierr);

        // check if everything has been set up correctly
        ierr=Assert(this->m_Optimizer!=NULL,"optimizer is not setup"); CHKERRQ(ierr);
        ierr=Assert(this->m_RegProblem!=NULL,"registration problem is not set up"); CHKERRQ(ierr);

        // set up synthetic test problem
        ierr=this->m_RegProblem->SetupSyntheticProb(this->m_ReferenceImage,this->m_TemplateImage); CHKERRQ(ierr);

    }

    // make sure images have not been set
    ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

    // allocate multilevel pyramid for reference image
    if (this->m_Opt->GetVerbosity() > 1){
        ierr=DbgMsg("setting up reference image multilevel pyramid"); CHKERRQ(ierr);
    }
    if(this->m_ReferencePyramid==NULL){
        try{this->m_ReferencePyramid = new MultiLevelPyramid(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // setup multilevel pyramid for reference image
    ierr=this->m_ReferencePyramid->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
    ierr=this->m_ReferencePyramid->DoSetup(this->m_ReferenceImage); CHKERRQ(ierr);

    // allocate multilevel pyramid for template image
    if (this->m_Opt->GetVerbosity() > 1){
        ierr=DbgMsg("setting up template image multilevel pyramid"); CHKERRQ(ierr);
    }
    if(this->m_TemplatePyramid==NULL){
        try{this->m_TemplatePyramid = new MultiLevelPyramid(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    // setup multilevel pyramid for template image
    ierr=this->m_TemplatePyramid->SetPreProc(this->m_PreProc); CHKERRQ(ierr);
    ierr=this->m_TemplatePyramid->DoSetup(this->m_TemplateImage); CHKERRQ(ierr);

    // get grid size
    for (int i = 0; i < 3; ++i){
        nx[i] = this->m_Opt->GetGridContPara().nx[0][i];
    }
    ierr=this->m_Opt->GetSizes(nx,nl,ng); CHKERRQ(ierr);

    // TODO: allow for warm start
    try{v = new VecField(nl,ng);}
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr=v->SetValue(0.0); CHKERRQ(ierr);


    // reset tolerance for gradient (optimization); we do not want
    // to solve as accurately when we solve on the coarse grid
    greltol = this->m_Opt->GetOptPara().tol[2];

    if (greltol < 1E-2){

        if (this->m_Opt->GetVerbosity() > 1){
            ss  << std::scientific
                << "increasing tolerance for gradient: "
                << greltol << " >> " << tolscale*greltol;
            ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();
        }

        this->m_Opt->SetOptTol(2,tolscale*greltol);

    }

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // get number of levels
    nlevels = this->m_Opt->GetGridContPara().nlevels;

    // run multi-level solver
    level=0;
    while (level < nlevels){

        // get number of grid points for current level
        for (int i = 0; i < 3; ++i){
            nx[i] = this->m_Opt->GetGridContPara().nx[level][i];
        }
        nl = this->m_Opt->GetGridContPara().nlocal[level];
        ng = this->m_Opt->GetGridContPara().nglobal[level];

        // display user message
        ss << std::scientific << "level " << std::setw(3) << level
           <<"    nx=("<< nx[0]<<","<< nx[1]<<","<< nx[2]
           << "); (nl,ng)=("<< nl << "," << ng << ")";
         ierr=this->DispLevelMsg(ss.str(),rank); CHKERRQ(ierr);
         ss.str( std::string() ); ss.clear();

        // get the individual images from the pyramid
        ierr=this->m_ReferencePyramid->GetLevel(&mR,level); CHKERRQ(ierr);
        ierr=this->m_TemplatePyramid->GetLevel(&mT,level); CHKERRQ(ierr);

        // initialize
        for (int i=0; i<3; ++i){
            this->m_Opt->SetNumGridPoints(i,nx[i]);
        }
        ierr=this->m_Opt->DoSetup(false); CHKERRQ(ierr);

        // store intermediate results
        if (this->m_Opt->GetReadWriteFlags().iterates){

            ss << "reference-image-level=" << level << ".nii.gz";
            ierr=this->m_ReadWrite->Write(mR,ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();

            ss << "template-image-level=" << level << ".nii.gz";
            ierr=this->m_ReadWrite->Write(mT,ss.str()); CHKERRQ(ierr);
            ss.str( std::string() ); ss.clear();

        }

        // do the setup
        ierr=this->SetupSolver(); CHKERRQ(ierr);

        // set images
        ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

        // compute initial gradient, objective and
        // distance mesure for zero velocity field
        ierr=this->m_RegProblem->InitializeOptimization(); CHKERRQ(ierr);

        // set initial guess and registraiton problem
        ierr=this->m_Optimizer->SetInitialGuess(v); CHKERRQ(ierr);
        ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

        // reset tolerances
        if ( (level == (nlevels-1)) && (greltol < 1E-2) ){

            if (this->m_Opt->GetVerbosity() > 1){
                ss  << std::scientific
                    << "reseting tolerance for gradient: "
                    << tolscale*greltol << " >> " << greltol;
                ierr=DbgMsg(ss.str()); CHKERRQ(ierr);
                ss.str( std::string() ); ss.clear();
            }

            this->m_Opt->SetOptTol(2,greltol);

        }

        // run the optimizer
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        // get and parse solution
        ierr=this->m_Optimizer->GetSolution(xstar); CHKERRQ(ierr);
        ierr=v->SetComponents(xstar); CHKERRQ(ierr);

        ++level; // increment iterator

        if (level < nlevels){

            ierr=this->ProlongVelocityField(v,level); CHKERRQ(ierr);

            if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }
            if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }

        }

    }

    // get the solution
    ierr=Assert(this->m_Solution!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=this->m_Solution->Copy(v); CHKERRQ(ierr);

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    if (v!=NULL){ delete v; v = NULL;};
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); mR=NULL; }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); mT=NULL; }

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief prolong velocity field
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ProlongVelocityField"
PetscErrorCode RegistrationInterface::ProlongVelocityField(VecField*& v, int level)
{
    PetscErrorCode ierr;
    IntType nx_f[3],nx_c[3],nl_f,ng_f;
    VecField *v_f=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=Assert(v!=NULL,"null pointer"); CHKERRQ(ierr);

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("prolonging velocity field"); CHKERRQ(ierr);
    }

    // set up preprocessing
    if(this->m_PreProc==NULL){
        try{this->m_PreProc = new PreProcReg(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    this->m_PreProc->ResetGridChangeOps(true);

    // get number of grid points for current level
    for (int i = 0; i<3; ++i){
        nx_f[i] = this->m_Opt->GetGridContPara().nx[level  ][i];
        nx_c[i] = this->m_Opt->GetGridContPara().nx[level-1][i];
    }

    // get number of points to allocate
    ierr=this->m_Opt->GetSizes(nx_f,nl_f,ng_f); CHKERRQ(ierr);

    // allocate container for velocity field
    try{ v_f = new reg::VecField(nl_f,ng_f); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }

    // apply prolongation operator
    ierr=this->m_PreProc->Prolong(v_f,v,nx_f,nx_c); CHKERRQ(ierr);

    // allocate container for velocity field
    delete v; v = NULL;
    try{ v = new reg::VecField(nl_f,ng_f); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
    }
    ierr=v->Copy(v_f); CHKERRQ(ierr);

    if (v_f!=NULL) { delete v_f; v_f=NULL; }

    this->m_PreProc->ResetGridChangeOps(false);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief finalize optimization (displays information for user)
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Finalize"
PetscErrorCode RegistrationInterface::Finalize()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    // finalize optimizer (show tao output)
    ierr=this->m_Optimizer->Finalize(); CHKERRQ(ierr);

    // display time to solution
    ierr=this->m_Opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief run postprocessing of input data
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "RunPostProcessing"
PetscErrorCode RegistrationInterface::RunPostProcessing()
{
    PetscErrorCode ierr;
    Vec mR=NULL,mT=NULL;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=this->SetupRegProblem(); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!=NULL,"null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr=Assert(this->m_Solution!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_TemplateImage!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(this->m_ReferenceImage!=NULL,"null pointer"); CHKERRQ(ierr);

    // allocate image containers
    ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);
    ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);

    // allocate preprocessing class
    if(this->m_PreProc==NULL){
        try{this->m_PreProc = new PreProcReg(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    if (this->m_Opt->GetRegFlags().smoothingenabled){
        // apply smoothing
        ierr=this->m_PreProc->ApplySmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_PreProc->ApplySmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);
    }
    else{
        // copy input images
        ierr=VecCopy(this->m_ReferenceImage,mR); CHKERRQ(ierr);
        ierr=VecCopy(this->m_TemplateImage,mT); CHKERRQ(ierr);
    }

    // set reference and template images
    ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
    ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    // compute stuff
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vectors
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}



/********************************************************************
 * @brief compute deformation map or deformation gradient
 ********************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ComputeDefFields"
PetscErrorCode RegistrationInterface::ComputeDefFields()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    this->m_Opt->Enter(__FUNCT__);

    ierr=this->SetupRegProblem(); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!=NULL,"null pointer"); CHKERRQ(ierr);

    // user needs to set template and reference image and the solution
    ierr=Assert(this->m_Solution!=NULL,"null pointer"); CHKERRQ(ierr);

    // compute stuff
    ierr=this->m_RegProblem->SetControlVariable(this->m_Solution); CHKERRQ(ierr);

    if (this->m_Opt->GetReadWriteFlags().detdefgrad){
        ierr=Msg("computing determinant of deformation gradient"); CHKERRQ(ierr);
        ierr=this->m_RegProblem->ComputeDetDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().defgrad){
        ierr=Msg("computing deformation gradient"); CHKERRQ(ierr);
        ierr=this->m_RegProblem->ComputeDefGrad(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().defmap){
        ierr=Msg("computing deformation map"); CHKERRQ(ierr);
        ierr=this->m_RegProblem->ComputeDeformationMap(true); CHKERRQ(ierr);
    }
    if (this->m_Opt->GetReadWriteFlags().deffield){
        ierr=Msg("computing displacement field"); CHKERRQ(ierr);
        ierr=this->m_RegProblem->ComputeDisplacementField(true); CHKERRQ(ierr);
    }

    this->m_Opt->Exit(__FUNCT__);

    PetscFunctionReturn(0);

}




} // end of name space




#endif // _REGISTRATIONINTERFACE_CPP_
