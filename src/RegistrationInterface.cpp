
/**
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
    this->m_IO=NULL;
    this->m_Optimizer=NULL;
    this->m_RegProblem=NULL;
    this->m_Prepoc = NULL;

    this->m_Solution=NULL;
    this->m_InitialGuess=NULL;
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
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (this->m_RegProblem != NULL){
        delete this->m_RegProblem;
        this->m_RegProblem = NULL;
    }
    if (this->m_Optimizer != NULL){
        delete this->m_Optimizer;
        this->m_Optimizer = NULL;
    }
    if (this->m_Prepoc != NULL){
        delete this->m_Prepoc;
        this->m_Prepoc = NULL;
    }

    if (this->m_InitialGuess != NULL){
        ierr=VecDestroy(&this->m_InitialGuess); CHKERRQ(ierr);
        this->m_InitialGuess=NULL;
    }
    if (this->m_Solution != NULL){
        ierr=VecDestroy(&this->m_Solution); CHKERRQ(ierr);
        this->m_Solution=NULL;
    }

    PetscFunctionReturn(0);
}



/********************************************************************
 * @brief set initial guess
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetInitialGuess"
PetscErrorCode RegistrationInterface::SetInitialGuess(Vec x0)
{
    PetscErrorCode ierr;
    IntType nlu,ngu;

    PetscFunctionBegin;

    // the input better is not zero
    ierr=Assert(x0!=NULL,"null pointer"); CHKERRQ(ierr);

    // compute number of unknowns
    nlu = 3*this->m_Opt->GetNLocal();
    ngu = 3*this->m_Opt->GetNGlobal();

    // if we have not setup initial guess, do so
    if (this->m_InitialGuess == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_InitialGuess,nlu,ngu); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_InitialGuess); CHKERRQ(ierr);
    }

    ierr=VecCopy(x0,this->m_InitialGuess); CHKERRQ(ierr);


    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set read write operator
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetIO"
PetscErrorCode RegistrationInterface::SetIO(ReadWriteReg* io)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(io != NULL, "null pointer"); CHKERRQ(ierr);
    this->m_IO = io;

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

    ierr=Assert(mR != NULL, "input reference image is null pointer"); CHKERRQ(ierr);

    this->m_ReferenceImage = mR;
    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set template image (i.e., the template image)
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReferenceImage"
PetscErrorCode RegistrationInterface::SetTemplateImage(Vec mT)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(mT != NULL, "input reference image is null pointer"); CHKERRQ(ierr);

    this->m_TemplateImage = mT;

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
    IntType nlu,ngu;
    PetscFunctionBegin;

    // reset registration problem
    if (this->m_RegProblem != NULL){
        delete this->m_RegProblem; this->m_RegProblem = NULL;
    }
    // reset optimization problem
    if (this->m_Optimizer != NULL){
        delete this->m_Optimizer; this->m_Optimizer = NULL;
    }

    // allocate class for io
    try{ this->m_Optimizer = new OptimizerType(this->m_Opt); }
    catch (std::bad_alloc&){
        ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
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
        try{ this->m_RegProblem = new OptimalControlRegistrationIC(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }
    else{
        try{ this->m_RegProblem = new OptimalControlRegistrationRelaxedIC(this->m_Opt); }
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    ierr=Assert(this->m_IO!= NULL, "io is null"); CHKERRQ(ierr);
    ierr=this->m_RegProblem->SetIO(this->m_IO); CHKERRQ(ierr);

    nlu = 3*this->m_Opt->GetNLocal();
    ngu = 3*this->m_Opt->GetNGlobal();

    if (this->m_Solution != NULL){
        ierr=VecDestroy(&this->m_Solution); CHKERRQ(ierr);
        this->m_Solution = NULL;
    }
    ierr=VecCreate(PETSC_COMM_WORLD,&this->m_Solution); CHKERRQ(ierr);
    ierr=VecSetSizes(this->m_Solution,nlu,ngu); CHKERRQ(ierr);
    ierr=VecSetFromOptions(this->m_Solution); CHKERRQ(ierr);

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
    int rank;
    IntType nlu,ngu;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // compute number of unknowns
    nlu = 3*this->m_Opt->GetNLocal();
    ngu = 3*this->m_Opt->GetNGlobal();

    // create all zero initial guess, if it has not been set already
    if (this->m_InitialGuess == NULL){
        ierr=VecCreate(PETSC_COMM_WORLD,&this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSetSizes(this->m_InitialGuess,nlu,ngu); CHKERRQ(ierr);
        ierr=VecSetFromOptions(this->m_InitialGuess); CHKERRQ(ierr);
        ierr=VecSet(this->m_InitialGuess,0.0); CHKERRQ(ierr);
    }

    ierr=Msg("starting optimization"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=PetscPrintf(PETSC_COMM_WORLD," %s  %-20s %-20s %-20s %-20s %-20s\n",
                                    "iter","objective (rel)","mismatch (rel)",
                                    "||gradient||_2,rel","||gradient||_2","step"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;

    // switch between all the options we have
    // for solving the optimization probelm
    if( this->m_Opt->DoRegParaCont() ){

        // run different flavours of parameter continuation
        ierr=this->RunSolverRegParaCont(); CHKERRQ(ierr);

    }
    else if( this->m_Opt->DoScaleCont() ){

        // run different flavours of parameter continuation
        ierr=this->RunSolverScaleCont(); CHKERRQ(ierr);

    }
    else{ ierr=this->RunSolver(); CHKERRQ(ierr); }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg("optimization done"); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;

    // finalize
    ierr=this->m_Optimizer->Finalize(); CHKERRQ(ierr);
    ierr=this->m_Opt->DisplayTimeToSolution(); CHKERRQ(ierr);


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

    // do the setup
    ierr=this->SetupRegProblem(); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!= NULL, "registration problem is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Optimizer!= NULL, "optimizer is null"); CHKERRQ(ierr);

    // presmoothing, if necessary
    if (this->m_Opt->ReadImagesFromFile()){

        ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
        ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

        // allocate
        ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);
        ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);

        if(this->m_Prepoc==NULL){
            try{this->m_Prepoc = new PreProcessingRegistration(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

        // rescale images
        ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
        ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

        ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{ ierr=this->m_RegProblem->SetupSyntheticProb(); CHKERRQ(ierr); }

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // init solver
    ierr=this->m_Optimizer->SetInitialGuess(this->m_InitialGuess); CHKERRQ(ierr);
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // run the optimization
    ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vector
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

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
    Vec mT=NULL,mR=NULL,x=NULL;
    PetscFunctionBegin;

    // do the setup
    ierr=this->SetupRegProblem(); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!= NULL, "registration problem is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Optimizer!= NULL, "optimizer is null"); CHKERRQ(ierr);

    // presmoothing, if necessary
    if (this->m_Opt->ReadImagesFromFile()){

        ierr=Assert(this->m_TemplateImage!=NULL,"template image is null"); CHKERRQ(ierr);
        ierr=Assert(this->m_ReferenceImage!=NULL,"reference image is null"); CHKERRQ(ierr);

        // allocate
        ierr=VecDuplicate(this->m_TemplateImage,&mT); CHKERRQ(ierr);
        ierr=VecDuplicate(this->m_ReferenceImage,&mR); CHKERRQ(ierr);

        if(this->m_Prepoc==NULL){
            try{this->m_Prepoc = new PreProcessingRegistration(this->m_Opt);}
            catch (std::bad_alloc&){
                ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
            }
        }

        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

        // rescale images
        ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
        ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

        ierr=this->m_RegProblem->SetReferenceImage(mR); CHKERRQ(ierr);
        ierr=this->m_RegProblem->SetTemplateImage(mT); CHKERRQ(ierr);

    }
    else{ ierr=this->m_RegProblem->SetupSyntheticProb(); CHKERRQ(ierr); }

    // reset all the clocks we have used so far
    ierr=this->m_Opt->ResetTimers(); CHKERRQ(ierr);
    ierr=this->m_Opt->ResetCounters(); CHKERRQ(ierr);

    // switch between the different strategies for
    // doing the parameter continuation (default one
    // is binary search)
    switch (this->m_Opt->GetRegParaContStrategy()){
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
#define __FUNCT__ "RunSolverRegParaBinarySearch"
PetscErrorCode RegistrationInterface::RunSolverRegParaContBinarySearch()
{
    PetscErrorCode ierr;
    int maxsteps,level,rank;
    bool stop,boundreached;
    std::stringstream ss;
    ScalarType beta,betamin,betascale,dbetascale,
                betastar,betahat,dbeta,dbetamin;
    Vec x;

    PetscFunctionBegin;

    ierr=Assert(this->m_Optimizer!=NULL,"optimizer is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!=NULL,"registration problem is null"); CHKERRQ(ierr);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    maxsteps  = this->m_Opt->GetMaxStepsParaCont();
    betamin   = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetBetaScaleParaCont();

    ierr=Assert(betascale < 1.0 && betascale > 0.0,"scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0 && betamin < 1.0,"lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess
    ierr=VecCopy(this->m_InitialGuess,this->m_Solution); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // initialize parameters
    beta = 1.0;
    betastar = beta;

    // reduce regularization parameter by one order of magnitude until
    // we hit tolerance
    stop=false; level = 0;
    while(level < maxsteps){

        this->m_Opt->SetRegularizationWeight(beta);

        ss << std::scientific << std::setw(3)
            << level <<" ( betav="<<beta
            <<"; betav*="<<betastar<<" )";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // set initial guess for current level
        ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        // get the solution
        ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);

        // check bounds on jacobian
        ierr=this->m_RegProblem->CheckBounds(x,stop); CHKERRQ(ierr);

        if (stop) break; // if bound reached go home

        // remember regularization parameter
        betastar = beta;

        // if we got here, the solution is valid
        ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

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
    dbetascale = this->m_Opt->GetDeltaBetaScaleParaCont();
    ierr=Assert(dbetascale<1.0 && dbetascale>0.0,"scale for delta betav not in (0,1)"); CHKERRQ(ierr);
    //update beta
    dbetamin = dbetascale*betastar;
    betahat = betascale*betastar;
    dbeta = (betastar-betahat)/2.0;
    beta  = betastar-dbeta;
    ++level;

    while(!stop){

        // set regularization parameter
        this->m_Opt->SetRegularizationWeight(beta);

        // display regularization parameter to user
        ss << std::setw(3) << level <<" ( betav="<<beta<<"; betav*="<<betastar<<" )";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // set initial guess for current level
        ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

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
            ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

        }

        // increase or reduce beta
        dbeta = (betastar - betahat)/2.0;
        beta  = betastar - dbeta;
        if (fabs(dbeta) < dbetamin){ stop = true; }

        ++level;

    }

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss <<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

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

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get parameters
    maxsteps  = this->m_Opt->GetMaxStepsParaCont();
    betamin   = this->m_Opt->GetBetaMinParaCont();
    betascale = this->m_Opt->GetBetaScaleParaCont();

    ierr=Assert(betascale < 1.0 && betascale > 0.0,"scale for beta not in (0,1)"); CHKERRQ(ierr);
    ierr=Assert(betamin > 0.0 && betamin < 1.0,"lower bound for beta in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess
    ierr=VecCopy(this->m_InitialGuess,this->m_Solution); CHKERRQ(ierr);

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
        this->m_Opt->SetRegularizationWeight(beta);

        // display message to user
        ss << std::scientific << std::setw(3) << level <<" (beta="<<beta<<")";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
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
        ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

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

    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss <<std::scientific<<"estimated regularization parameter betav="<<betastar;
    ierr=Msg(ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

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

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // get target regularization weight
    betastar=this->m_Opt->GetRegularizationWeight(0);
    ierr=Assert(betastar>0.0 && betastar<1.0,"target beta not in (0,1)"); CHKERRQ(ierr);

    // set optimization problem
    ierr=this->m_Optimizer->SetProblem(this->m_RegProblem); CHKERRQ(ierr);

    // set initial guess
    ierr=VecCopy(this->m_InitialGuess,this->m_Solution); CHKERRQ(ierr);

    // set initial guess for current level
    ierr=this->m_Optimizer->SetInitialGuess(this->m_Solution); CHKERRQ(ierr);

    // reduce regularization parameter
    level = 0; beta=1.0;
    while(beta > betastar){

        // set regularization weight
        this->m_Opt->SetRegularizationWeight(beta);

        // display message to user
        ss << std::scientific << std::setw(3) << level <<" (beta="<<beta<<"; beta*="<<betastar<<")";
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        beta /= static_cast<ScalarType>(10); // reduce beta
        ++level;

    } // parameter reduction

    beta = betastar;

    // set regularization weight
    this->m_Opt->SetRegularizationWeight(beta);

    // display message to user
    ss << std::scientific << std::setw(3) << level <<" (beta="<<beta<<"; beta*="<<betastar<<")";
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
    if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
    ss.str( std::string() ); ss.clear();

    // solve optimization problem for user defined regularization parameter
    ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

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
    ScalarType sigma[3];

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // set up preprocessing
    if(this->m_Prepoc==NULL){
        try{this->m_Prepoc = new PreProcessingRegistration(this->m_Opt);}
        catch (std::bad_alloc&){
            ierr=reg::ThrowError("allocation failed"); CHKERRQ(ierr);
        }
    }

    // do the setup
    ierr=this->SetupRegProblem(); CHKERRQ(ierr);
    ierr=Assert(this->m_RegProblem!= NULL,"registration problem is null"); CHKERRQ(ierr);
    ierr=Assert(this->m_Optimizer!= NULL,"optimizer is null"); CHKERRQ(ierr);

    // set up synthetic problem if we did not read images
    if (!this->m_Opt->ReadImagesFromFile()){
        // TODO: handle synthetic test problem
        ierr=this->m_RegProblem->SetupSyntheticProb(); CHKERRQ(ierr);
    }

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
    ierr=this->m_Optimizer->SetInitialGuess(this->m_InitialGuess); CHKERRQ(ierr);

    level=0;
    while (level < maxlevel){

        // this resets the optimizer to assume that we're in first iteration
        // every time we solve the problem; TODO: add initialization stage
        // for zero velocity field (compute gradient) to refer warm starts
        // and every thing else to the initial gradient
        this->m_RegProblem->InFirstIteration(true);

        for (int i=0; i < 3; ++i){
            sigma[i] = this->m_Opt->GetScaleContSigma(i,level);
            this->m_Opt->SetSigma(i,sigma[i]);
        }

        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mR,this->m_ReferenceImage); CHKERRQ(ierr);
        ierr=this->m_Prepoc->ApplyGaussianSmoothing(mT,this->m_TemplateImage); CHKERRQ(ierr);

        // rescale images
        ierr=Rescale(mR,0.0,1.0); CHKERRQ(ierr);
        ierr=Rescale(mT,0.0,1.0); CHKERRQ(ierr);

        // display message to user
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss << std::scientific << std::setw(3) << level <<" sigma=("
            <<sigma[0]<<","<<sigma[1]<<","<<sigma[2]<<")";
        ierr=Msg("level " + ss.str()); CHKERRQ(ierr);
        if (rank == 0) std::cout << std::string(this->m_Opt->GetLineLength(),'-') << std::endl;
        ss.str( std::string() ); ss.clear();

        // run the optimization
        ierr=this->m_Optimizer->Run(); CHKERRQ(ierr);

        ++level;

    }

    // get the solution
    ierr=this->m_Optimizer->GetSolution(x); CHKERRQ(ierr);
    ierr=VecCopy(x,this->m_Solution); CHKERRQ(ierr);

    // wrap up
    ierr=this->m_RegProblem->Finalize(this->m_Solution); CHKERRQ(ierr);

    // destroy vector
    if (mR!=NULL){ ierr=VecDestroy(&mR); CHKERRQ(ierr); }
    if (mT!=NULL){ ierr=VecDestroy(&mT); CHKERRQ(ierr); }

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

    // display time to solution
    ierr=this->m_Opt->DisplayTimeToSolution(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}




} // end of name space




#endif // _REGISTRATIONINTERFACE_CPP_
