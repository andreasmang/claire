#ifndef _PREPROCREG_CPP_
#define _PREPROCREG_CPP_


#include "PreProcReg.hpp"


namespace reg
{




/********************************************************************
 * @brief default constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcReg"
PreProcReg::PreProcReg()
{
    this->Initialize();
}




/********************************************************************
 * @brief constructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "PreProcReg"
PreProcReg::PreProcReg(RegOpt* opt)
{
    this->Initialize();
    this->m_Opt = opt;
}




/********************************************************************
 * @brief default deconstructor
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "~PreProcReg"
PreProcReg::~PreProcReg()
{
    this->ClearMemory();
}




/********************************************************************
 * @brief initialize
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Initialize"
PetscErrorCode PreProcReg::Initialize()
{
    this->m_Opt = NULL;

    this->m_ResetGridChangeOperators=false;

    this->m_ReadWrite = NULL;

    this->m_XHatFine = NULL;
    this->m_XHatCoarse = NULL;
    this->m_FFTFinePlan = NULL;
    this->m_FFTCoarsePlan = NULL;

    this->m_xhat = NULL;
    this->m_yhat = NULL;

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief clear memory
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ClearMemory"
PetscErrorCode PreProcReg::ClearMemory()
{
    //PetscErrorCode ierr;
    if(this->m_xhat!=NULL){
        accfft_free(this->m_xhat);
        this->m_xhat = NULL;
    }
    if(this->m_yhat!=NULL){
        accfft_free(this->m_yhat);
        this->m_yhat = NULL;
    }

    if(this->m_XHatFine != NULL){
        accfft_free(this->m_XHatFine);
        this->m_XHatFine=NULL;
    }
    if(this->m_XHatCoarse != NULL){
        accfft_free(this->m_XHatCoarse);
        this->m_XHatCoarse=NULL;
    }

    if(this->m_FFTFinePlan != NULL){
        accfft_destroy_plan(this->m_FFTFinePlan);
        this->m_FFTFinePlan=NULL;
    }
    if(this->m_FFTCoarsePlan != NULL){
        accfft_destroy_plan(this->m_FFTCoarsePlan);
        this->m_FFTCoarsePlan=NULL;
    }

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief set read/write object for data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetReadWrite"
PetscErrorCode PreProcReg::SetReadWrite(PreProcReg::ReadWriteType* readwrite)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr=Assert(readwrite!= NULL,"null pointer"); CHKERRQ(ierr);
    this->m_ReadWrite=readwrite;

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief do setup for applying prolongation and restriction
 * operators
 * @param nx_f grid size on fine grid
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupGridChangeOperators"
PetscErrorCode PreProcReg::SetupGridChangeOperators(IntType* nx_f,IntType* nx_c)
{
    PetscErrorCode ierr;
    IntType nalloc_c,nalloc_f;
    int _nx_f[3],_ostart_f[3],_osize_f[3],_isize_f[3],_istart_f[3],
        _nx_c[3],_ostart_c[3],_osize_c[3],_isize_c[3],_istart_c[3];
    ScalarType *p_xfd=NULL,*p_xcd=NULL;
    Complex *p_xfdhat=NULL,*p_xcdhat=NULL;
    MPI_Comm mpicomm;
    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("setting up grid change operators"); CHKERRQ(ierr);
    }

    // parse input sizes
    _nx_f[0] = static_cast<int>(nx_f[0]);
    _nx_f[1] = static_cast<int>(nx_f[1]);
    _nx_f[2] = static_cast<int>(nx_f[2]);

    _nx_c[0] = static_cast<int>(nx_c[0]);
    _nx_c[1] = static_cast<int>(nx_c[1]);
    _nx_c[2] = static_cast<int>(nx_c[2]);

    this->m_FFTFineScale = 1.0;
    this->m_FFTCoarseScale = 1.0;
    for (int i=0; i < 3; ++i){
        this->m_FFTFineScale *= static_cast<ScalarType>(nx_f[i]);
        this->m_FFTCoarseScale *= static_cast<ScalarType>(nx_c[i]);
    }
    this->m_FFTFineScale = 1.0/this->m_FFTFineScale;
    this->m_FFTCoarseScale = 1.0/this->m_FFTCoarseScale;


    // get communicator
    mpicomm = this->m_Opt->GetFFT().mpicomm;

    if ( this->m_ResetGridChangeOperators ){

        if (this->m_XHatCoarse != NULL){
            std::cout<< "deleting x hat coarse"<<std::endl;
            accfft_free(this->m_XHatCoarse);
            this->m_XHatCoarse=NULL;
        }
        if (this->m_XHatFine != NULL){
            std::cout<< "deleting x hat fine"<<std::endl;
            accfft_free(this->m_XHatFine);
            this->m_XHatFine=NULL;
        }

        if(this->m_FFTFinePlan != NULL){
            std::cout<< "deleting fft plan fine"<<std::endl;
            accfft_destroy_plan(this->m_FFTFinePlan);
            this->m_FFTFinePlan=NULL;
            accfft_cleanup();
        }

        if(this->m_FFTCoarsePlan != NULL){
            std::cout<< "deleting fft plan coarse"<<std::endl;
            accfft_destroy_plan(this->m_FFTCoarsePlan);
            this->m_FFTCoarsePlan=NULL;
            accfft_cleanup();
        }
    }

    nalloc_c = accfft_local_size_dft_r2c(_nx_c,_isize_c,_istart_c,_osize_c,_ostart_c,mpicomm);
    nalloc_f = accfft_local_size_dft_r2c(_nx_f,_isize_f,_istart_f,_osize_f,_ostart_f,mpicomm);

    for (int i = 0; i < 3; ++i){
        this->m_osize_c[i] = static_cast<IntType>(_osize_c[i]);
        this->m_osize_f[i] = static_cast<IntType>(_osize_f[i]);

        this->m_ostart_c[i] = static_cast<IntType>(_ostart_c[i]);
        this->m_ostart_f[i] = static_cast<IntType>(_ostart_f[i]);
    }

    if (this->m_XHatCoarse == NULL){
        std::cout<< "allocating x hat coarse"<<std::endl;
        this->m_XHatCoarse = (ScalarTypeFD*)accfft_alloc(nalloc_c);
    }

    if (this->m_XHatFine == NULL){
        std::cout<< "allocating x hat fine"<<std::endl;
        this->m_XHatFine = (ScalarTypeFD*)accfft_alloc(nalloc_f);
    }

    // allocate plan for fine grid
    if (this->m_FFTFinePlan == NULL){

        p_xfd = (ScalarType*)accfft_alloc(nalloc_f);
        p_xfdhat = (Complex*)accfft_alloc(nalloc_f);

        std::cout<< "allocating fft plan fine"<<std::endl;
        this->m_FFTFinePlan = accfft_plan_dft_3d_r2c(_nx_f,p_xfd,(double*)p_xfdhat,mpicomm,ACCFFT_MEASURE);

        accfft_free(p_xfd); p_xfd=NULL;
        accfft_free(p_xfdhat); p_xfdhat=NULL;
    }

    // allocate plan for coarse grid
    if (this->m_FFTCoarsePlan == NULL){

        p_xcd = (ScalarType*)accfft_alloc(nalloc_c);
        p_xcdhat = (Complex*)accfft_alloc(nalloc_c);

        std::cout<< "allocating fft plan coarse"<<std::endl;
        this->m_FFTCoarsePlan = accfft_plan_dft_3d_r2c(_nx_c,p_xcd,(double*)p_xcdhat,mpicomm,ACCFFT_MEASURE);

        accfft_free(p_xcd); p_xcd=NULL;
        accfft_free(p_xcdhat); p_xcdhat=NULL;
    }

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("setup of grid change operators done"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief restrict vector field
 * @param v input vector field
 * @param vcoarse output vector field v_c = R[v]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(VecField* vcoarse, VecField* vfine, IntType* nx_c, IntType* nx_f)
{
    PetscErrorCode ierr;

    ierr=Assert(vfine!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(vcoarse!=NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Restrict(&vcoarse->m_X1,vfine->m_X1,nx_c,nx_f); CHKERRQ(ierr);
    ierr=this->Restrict(&vcoarse->m_X2,vfine->m_X2,nx_c,nx_f,false); CHKERRQ(ierr);
    ierr=this->Restrict(&vcoarse->m_X3,vfine->m_X3,nx_c,nx_f,false); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief restrict data
 * @param x input vector
 * @param xcoarse output vector xcoarse = R[x]
 * @param nx_c number of grid points on coarse grid
 * @param dosetup flag to identify if we have to do the setup step;
 * this is essentially to prevent an additional setup step if we
 * apply this function to each component of a vector field, or a
 * time dependend field; if the parameter is not set, it is true
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Restrict"
PetscErrorCode PreProcReg::Restrict(Vec* x_c, Vec x_f, IntType* nx_c, IntType* nx_f, bool dosetup)
{
    PetscErrorCode ierr;
    ScalarType *p_xf=NULL,*p_xc=NULL,scale;
    IntType n;
    int rank;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying restriction operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr=Assert(x_f!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(x_c!=NULL,"null pointer"); CHKERRQ(ierr);

    ierr=this->SetupGridChangeOperators(nx_f,nx_c); CHKERRQ(ierr);

    // set freqencies to zero
    for (IntType i1_c = 0; i1_c < this->m_osize_c[0]; ++i1_c){
        for (IntType i2_c = 0; i2_c < this->m_osize_c[1]; ++i2_c){
            for (IntType i3_c = 0; i3_c < this->m_osize_c[2]; ++i3_c){

                IntType i = GetLinearIndex(i1_c,i2_c,i3_c,this->m_osize_c);

                this->m_XHatCoarse[i][0] = 0.0;
                this->m_XHatCoarse[i][1] = 0.0;

            }
        }
    }

    // compute fft of data on fine grid
    ierr=VecGetArray(x_f,&p_xf); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_FFTFinePlan,p_xf,this->m_XHatFine,timer);
    ierr=VecRestoreArray(x_f,&p_xf); CHKERRQ(ierr);

    // get size
    if (dosetup) ierr=this->SetupRestriction(nx_f,nx_c); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTFineScale;

    // get number of entries we are going to assign
    n = this->m_IndicesC[rank].size();
    ierr=Assert(n==this->m_IndicesF[rank].size(),"size error"); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType lf,lc;
#pragma omp for
    for (IntType i = 0; i < n; ++i){

        lf = this->m_IndicesF[rank][i];
        lc = this->m_IndicesC[rank][i];

        this->m_XHatCoarse[lc][0] = scale*this->m_XHatFine[lf][0];
        this->m_XHatCoarse[lc][1] = scale*this->m_XHatFine[lf][1];

    }

} // pragma omp parallel


    ierr=VecGetArray(*x_c,&p_xc); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_FFTCoarsePlan,this->m_XHatCoarse,p_xc,timer);
    ierr=VecRestoreArray(*x_c,&p_xc); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(timer);


    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for applying restriction operator
 * @param nx_c grid size on coarse grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupRestriction"
PetscErrorCode PreProcReg::SetupRestriction(IntType* nx_f,IntType* nx_c)
{
    PetscErrorCode ierr;
    int rank,nprocs,nowned,nsend,nprocessed,xrank,c_grid[2],p1,p2;
    IntType oend_c[3],osizex2,osizex3,li_f,li_c,i1_f,i2_f,i3_f,i1_c,i2_c,i3_c;
    ScalarType k1_f,k2_f,k3_f,k1_c,k2_c,k3_c,nxhalf_c[3];
    bool owned;

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(nx_c[0] <= nx_f[0],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[1] <= nx_f[1],"grid size in restriction wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_c[2] <= nx_f[2],"grid size in restriction wrong"); CHKERRQ(ierr);

    // allocate if necessary
    if (this->m_IndicesC.empty()){
        this->m_IndicesC.resize(nprocs);
    }
    if (this->m_IndicesF.empty()){
        this->m_IndicesF.resize(nprocs);
    }

    for (int i = 0; i < nprocs; ++i){
        if (!this->m_IndicesC[i].empty()){
            this->m_IndicesC[i].clear();
        }
        if (!this->m_IndicesF[i].empty()){
            this->m_IndicesF[i].clear();
        }
    }

    for(int i = 0; i < 3; ++i){
        nxhalf_c[i] = std::floor(static_cast<ScalarType>(nx_c[i])/2.0);
        oend_c[i]   = this->m_ostart_c[i] + this->m_osize_c[i];
    }

    // get cartesian grid (MPI)
    c_grid[0] = this->m_Opt->GetNetworkDims(0);
    c_grid[1] = this->m_Opt->GetNetworkDims(1);

    osizex2 = std::ceil( static_cast<ScalarType>(nx_c[1])/static_cast<ScalarType>(c_grid[0]));
    osizex3 = std::ceil( (static_cast<ScalarType>(nx_c[2])/2.0 + 1.0)/static_cast<ScalarType>(c_grid[1]));

    // for all points on fine grid
    nowned=0; nsend=0; nprocessed=0;
    for (i1_f = 0; i1_f < this->m_osize_f[0]; ++i1_f){ // x1
        for (i2_f = 0; i2_f < this->m_osize_f[1]; ++i2_f){ // x2
            for (i3_f = 0; i3_f < this->m_osize_f[2]; ++i3_f){ // x3

                // compute wave number index on fine grid
                k1_f = static_cast<IntType>(i1_f + this->m_ostart_f[0]);
                k2_f = static_cast<IntType>(i2_f + this->m_ostart_f[1]);
                k3_f = static_cast<IntType>(i3_f + this->m_ostart_f[2]);

                // only if current fourier entry is represented in
                // spectral domain of coarse grid
                if (  ( k1_f <= nxhalf_c[0] || k1_f > (nx_f[0] - nxhalf_c[0]) )
                   && ( k2_f <= nxhalf_c[1] || k2_f > (nx_f[1] - nxhalf_c[1]) )
                   && ( k3_f <= nxhalf_c[2] || k3_f > (nx_f[2] - nxhalf_c[2]) ) ){

                    ++nprocessed;

                    // get wave number index on coarse grid
                    k1_c = k1_f <= nxhalf_c[0] ? k1_f : nx_c[0] - nx_f[0] + k1_f;
                    k2_c = k2_f <= nxhalf_c[1] ? k2_f : nx_c[1] - nx_f[1] + k2_f;
                    k3_c = k3_f <= nxhalf_c[2] ? k3_f : nx_c[2] - nx_f[2] + k3_f;

                    if ( (k1_c < 0.0) || (k2_c < 0.0) || (k3_c < 0.0) ){
                        std::cout<<"index out of bounds (smaller than zero)"<<std::endl;
                    }
                    if ( (k1_c > nx_c[0]) || (k2_c > nx_c[1]) || (k3_c > nx_c[2]) ){
                        std::cout<<"index out of bounds (larger than nx)"<<std::endl;
                    }

                    owned=true;
                    if ( (k1_c < this->m_ostart_c[0]) || (k1_c >= oend_c[0]) ) owned = false;
                    if ( (k2_c < this->m_ostart_c[1]) || (k2_c >= oend_c[1]) ) owned = false;
                    if ( (k3_c < this->m_ostart_c[2]) || (k3_c >= oend_c[2]) ) owned = false;

                    // compute processor id
                    p1=static_cast<int>(k2_c/osizex2);
                    p2=static_cast<int>(k3_c/osizex3);

                    xrank = p1*c_grid[1] + p2;

                    if ( owned ){

                        // compute local index
                        i1_c = static_cast<IntType>(k1_c) - this->m_ostart_c[0];
                        i2_c = static_cast<IntType>(k2_c) - this->m_ostart_c[1];
                        i3_c = static_cast<IntType>(k3_c) - this->m_ostart_c[2];

                        if ( (i1_c >= this->m_osize_c[0]) || (i2_c >= this->m_osize_c[1]) || (i3_c >= this->m_osize_c[2]) ){
                            std::cout<<"index out of bounds (larger than osize)"<<std::endl;
                        }
                        if ( (i1_c < 0.0) || (i2_c < 0.0) || (i3_c < 0.0) ){
                            std::cout<<"index out of bounds (negative)"<<std::endl;
                        }
                        // compute flat index
                        li_c = GetLinearIndex(i1_c,i2_c,i3_c,this->m_osize_c);
                        li_f = GetLinearIndex(i1_f,i2_f,i3_f,this->m_osize_f);

                        // assign computed indices to array (for given rank)
                        this->m_IndicesC[rank].push_back(li_c);
                        this->m_IndicesF[rank].push_back(li_f);

                        // check if woned is really owned
                        if ( rank != xrank ){
                            std::cout<<"rank not owned: " << rank << " " << xrank <<std::endl;
                        }
                        ++nowned;
                    }
                    else{

                        if ( rank == xrank ){
                            std::cout<<" rank owned: " << rank << " " << xrank <<std::endl;
                        }

                        this->m_IndicesC[xrank].push_back(0);
                        this->m_IndicesF[xrank].push_back(0);

                        ++nsend;

                    }
                }

            } // i1
        } // i2
    } // i3

    // TODO: send data that does not belong to current proc to
    // other procs
    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief prolong vector field
 * @param vcoarse input vector field
 * @param vfine output vector field vfine = P[vcoarse]
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(VecField* v_f, VecField* v_c, IntType* nx_f, IntType* nx_c)
{
    PetscErrorCode ierr;

    ierr=Assert(v_f!=NULL, "null pointer"); CHKERRQ(ierr);
    ierr=Assert(v_c!=NULL, "null pointer"); CHKERRQ(ierr);

    ierr=this->Prolong(&v_f->m_X1,v_c->m_X1,nx_f,nx_c); CHKERRQ(ierr);
    ierr=this->Prolong(&v_f->m_X2,v_c->m_X2,nx_f,nx_c,false); CHKERRQ(ierr);
    ierr=this->Prolong(&v_f->m_X3,v_c->m_X3,nx_f,nx_c,false); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief prolong scalar field
 * @param x input vector
 * @param xfine output vector xfine = P[x]
 * @param nx_f number of grid points on fine grid
 * @param dosetup flag to identify if we have to do the setup step;
 * this is essentially to prevent an additional setup step if we
 * apply this function to each component of a vector field, or a
 * time dependend field; if the parameter is not set, it is true
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "Prolong"
PetscErrorCode PreProcReg::Prolong(Vec* x_f, Vec x_c, IntType* nx_f, IntType* nx_c, bool dosetup)
{
    PetscErrorCode ierr;
    ScalarType *p_xf=NULL,*p_xc=NULL,scale;
    IntType n;
    int rank;
    double timer[5]={0,0,0,0,0};

    PetscFunctionBegin;

    if (this->m_Opt->GetVerbosity() > 2){
        ierr=DbgMsg("applying prolongation operator"); CHKERRQ(ierr);
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    ierr=Assert(x_c!=NULL,"null pointer"); CHKERRQ(ierr);
    ierr=Assert(x_f!=NULL,"null pointer"); CHKERRQ(ierr);

    // allocate operators
    if (dosetup) ierr=this->SetupGridChangeOperators(nx_f,nx_c); CHKERRQ(ierr);

    // set freqencies to zero
    for (IntType i1_f = 0; i1_f < this->m_osize_f[0]; ++i1_f){
        for (IntType i2_f = 0; i2_f < this->m_osize_f[1]; ++i2_f){
            for (IntType i3_f = 0; i3_f < this->m_osize_f[2]; ++i3_f){

                IntType i = GetLinearIndex(i1_f,i2_f,i3_f,this->m_osize_f);

                this->m_XHatFine[i][0] = 0.0;
                this->m_XHatFine[i][1] = 0.0;

            }
        }
    }

    // get size
    if (dosetup) ierr=this->SetupProlongation(nx_f,nx_c); CHKERRQ(ierr);

    // compute fft of data on fine grid
    ierr=VecGetArray(x_c,&p_xc); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_FFTCoarsePlan,p_xc,this->m_XHatCoarse,timer);
    ierr=VecRestoreArray(x_c,&p_xc); CHKERRQ(ierr);

    // get grid sizes/fft scales
    scale = this->m_FFTCoarseScale;

    // get number of entries we are going to assign
    n = this->m_IndicesC[rank].size();
    ierr=Assert(n==this->m_IndicesF[rank].size(),"size error"); CHKERRQ(ierr);

#pragma omp parallel
{
    IntType lc,lf;
#pragma omp for
    for (IntType i = 0; i < n; ++i){

        lc = this->m_IndicesC[rank][i];
        lf = this->m_IndicesF[rank][i];

        this->m_XHatFine[lf][0] = scale*this->m_XHatCoarse[lc][0];
        this->m_XHatFine[lf][1] = scale*this->m_XHatCoarse[lc][1];

    }

} // pragma omp parallel


    ierr=VecGetArray(*x_f,&p_xf); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_FFTFinePlan,this->m_XHatFine,p_xf,timer);
    ierr=VecRestoreArray(*x_f,&p_xf); CHKERRQ(ierr);

    // set fft timers
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(0);
}




/********************************************************************
 * @brief do setup for applying prolongation operator
 * @param nx_f grid size on fine grid
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "SetupProlongation"
PetscErrorCode PreProcReg::SetupProlongation(IntType* nx_f, IntType* nx_c)
{
    PetscErrorCode ierr;
    int rank,nprocs,nowned,nsend,xrank,c_grid[2],p1,p2;
    IntType oend_f[3],osizex2,osizex3,li_c,li_f,
            i1_c,i2_c,i3_c, i1_f,i2_f,i3_f;
    ScalarType k1_c,k2_c,k3_c, k1_f,k2_f,k3_f, nxhalf_c[3];

    PetscFunctionBegin;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&nprocs);

    ierr=Assert(nx_f[0] >= nx_c[0],"grid size in prolongation wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_f[1] >= nx_c[1],"grid size in prolongation wrong"); CHKERRQ(ierr);
    ierr=Assert(nx_f[2] >= nx_c[2],"grid size in prolongation wrong"); CHKERRQ(ierr);

    // allocate if necessary
    if (this->m_IndicesC.empty()){
        this->m_IndicesC.resize(nprocs);
    }
    if (this->m_IndicesF.empty()){
        this->m_IndicesF.resize(nprocs);
    }

    for (int i = 0; i < nprocs; ++i){
        if (!this->m_IndicesC[i].empty()){
            this->m_IndicesC[i].clear();
        }
        if (!this->m_IndicesF[i].empty()){
            this->m_IndicesF[i].clear();
        }
    }

    for(int i = 0; i < 3; ++i){
        oend_f[i] = this->m_ostart_f[i] + this->m_osize_f[i];
        nxhalf_c[i] = std::floor(static_cast<ScalarType>(nx_c[i])/2.0);
    }

    // get cartesian grid (MPI)
    c_grid[0] = this->m_Opt->GetNetworkDims(0);
    c_grid[1] = this->m_Opt->GetNetworkDims(1);

    osizex2 = std::ceil( static_cast<ScalarType>(nx_f[1])/static_cast<ScalarType>(c_grid[0]));
    osizex3 = std::ceil( (static_cast<ScalarType>(nx_f[2])/2.0 + 1.0)/static_cast<ScalarType>(c_grid[1]));

    // for all points on coarse grid
    nowned=0; nsend=0;
    for (i1_c = 0; i1_c < this->m_osize_c[0]; ++i1_c){ // x1
        for (i2_c = 0; i2_c < this->m_osize_c[1]; ++i2_c){ // x2
            for (i3_c = 0; i3_c < this->m_osize_c[2]; ++i3_c){ // x3

                // compute wave number index on fine grid
                k1_c = static_cast<IntType>(i1_c + this->m_ostart_c[0]);
                k2_c = static_cast<IntType>(i2_c + this->m_ostart_c[1]);
                k3_c = static_cast<IntType>(i3_c + this->m_ostart_c[2]);

                // get wave number index on coarse grid
                k1_f = k1_c <= nxhalf_c[0] ? k1_c : nx_f[0] - nx_c[0] + k1_c;
                k2_f = k2_c <= nxhalf_c[1] ? k2_c : nx_f[1] - nx_c[1] + k2_c;
                k3_f = k3_c <= nxhalf_c[2] ? k3_c : nx_f[2] - nx_c[2] + k3_c;

                if ( (k1_f < 0.0) || (k2_f < 0.0) || (k3_f < 0.0) ){
                    std::cout<<"index out of bounds (smaller than zero)"<<std::endl;
                }
                if ( (k1_f > nx_f[0]) || (k2_f > nx_f[1]) || (k3_f > nx_f[2]) ){
                    std::cout<<"index out of bounds (larger than nx)"<<std::endl;
                }

                // compute processor id
                p1=static_cast<int>(k2_f/osizex2);
                p2=static_cast<int>(k3_f/osizex3);
                xrank = p1*c_grid[1] + p2;

                // only if current fourier entry is represented in
                // spectral domain of coarse grid
                if (  ( k1_f >= this->m_ostart_f[0] && k1_f < oend_f[0] )
                   && ( k2_f >= this->m_ostart_f[1] && k2_f < oend_f[1] )
                   && ( k3_f >= this->m_ostart_f[2] && k3_f < oend_f[2] ) ){

                    // compute local index
                    i1_f = static_cast<IntType>(k1_f) - this->m_ostart_f[0];
                    i2_f = static_cast<IntType>(k2_f) - this->m_ostart_f[1];
                    i3_f = static_cast<IntType>(k3_f) - this->m_ostart_f[2];

                    if ( (i1_f >= this->m_osize_f[0]) || (i2_f >= this->m_osize_f[1]) || (i3_f >= this->m_osize_f[2]) ){
                        std::cout<<"index out of bounds (larger than osize)"<<std::endl;
                    }
                    if ( (i1_f < 0.0) || (i2_f < 0.0) || (i3_f < 0.0) ){
                        std::cout<<"index out of bounds (negative)"<<std::endl;
                    }
                    // compute flat index
                    li_f = GetLinearIndex(i1_f,i2_f,i3_f,this->m_osize_f);
                    li_c = GetLinearIndex(i1_c,i2_c,i3_c,this->m_osize_c);

                    // assign computed indices to array (for given rank)
                    this->m_IndicesF[rank].push_back(li_f);
                    this->m_IndicesC[rank].push_back(li_c);

                    ++nowned;

                    // check if woned is really owned
                    if ( rank != xrank ){
                        std::cout<<"rank not owned: " << rank <<" "<< xrank <<std::endl;
                    }
                }
                else{

                    if ( rank == xrank ){
                        std::cout<<" rank owned: "<< rank <<" "<< xrank <<std::endl;
                    }

                    this->m_IndicesF[xrank].push_back(0);
                    this->m_IndicesC[xrank].push_back(0);

                    ++nsend;

                }

            } // i1
        } // i2
    } // i3

    // TODO: send data that does not belong to current proc to
    // other procs
    PetscFunctionReturn(0);

}




/********************************************************************
 * @brief apply gaussian smoothing operator to input data
 *******************************************************************/
#undef __FUNCT__
#define __FUNCT__ "ApplySmoothing"
PetscErrorCode PreProcReg::ApplySmoothing(Vec xsmooth, Vec x)
{
    PetscErrorCode ierr;
    IntType osize[3],ostart[3],nalloc;
    ScalarType hx[3],nx[3],c[3],scale;
    ScalarType *p_x=NULL, *p_xsmooth=NULL;
    double timer[5]={0,0,0,0,0};
    PetscFunctionBegin;

    ierr=Assert(x != NULL,"null pointer"); CHKERRQ(ierr);

    // get local pencil size and allocation size
    nalloc=this->m_Opt->GetFFT().nalloc;

    if(this->m_xhat == NULL){
        this->m_xhat=(ScalarTypeFD*)accfft_alloc(nalloc);
    }

    // get parameters
    for (int i = 0; i < 3; ++i){

        osize[i] = this->m_Opt->GetFFT().osize[i];
        ostart[i] = this->m_Opt->GetFFT().ostart[i];

        hx[i] = this->m_Opt->GetDomainPara().hx[i];
        nx[i] = static_cast<ScalarType>(this->m_Opt->GetDomainPara().nx[i]);

        // sigma is provided by user in # of grid points
        c[i] = this->m_Opt->GetSigma(i)*hx[i];
        c[i] *= c[i];

    }

    scale = this->m_Opt->ComputeFFTScale();

    // compute fft
    ierr=VecGetArray(x,&p_x); CHKERRQ(ierr);
    accfft_execute_r2c_t<ScalarType,ScalarTypeFD>(this->m_Opt->GetFFT().plan,p_x,this->m_xhat,timer);
    ierr=VecRestoreArray(x,&p_x); CHKERRQ(ierr);


#pragma omp parallel
{
#pragma omp for
    for (IntType i1 = 0; i1 < osize[0]; ++i1){ // x1
        for (IntType i2 = 0; i2 < osize[1]; ++i2){ // x2
            for (IntType i3 = 0; i3 < osize[2]; ++i3){ // x3

                // compute coordinates (nodal grid)
                ScalarType k1 = static_cast<ScalarType>(i1 + ostart[0]);
                ScalarType k2 = static_cast<ScalarType>(i2 + ostart[1]);
                ScalarType k3 = static_cast<ScalarType>(i3 + ostart[2]);

                // check if grid index is larger or smaller then
                // half of the total grid size
                bool flagx1 = (k1 <= nx[0]*0.5);
                bool flagx2 = (k2 <= nx[1]*0.5);
                bool flagx3 = (k3 <= nx[2]*0.5);

                k1 = flagx1 ? k1 : -nx[0] + k1;
                k2 = flagx2 ? k2 : -nx[1] + k2;
                k3 = flagx3 ? k3 : -nx[2] + k3;

                ScalarType sik = 0.5*((k1*k1*c[0]) + (k2*k2*c[1]) + (k3*k3*c[2]));
                sik = exp(-sik);

                // compute linear / flat index
                IntType li = GetLinearIndex(i1,i2,i3,osize);

                this->m_xhat[li][0] *= scale*sik;
                this->m_xhat[li][1] *= scale*sik;

            } // i1
        } // i2
    } // i3

} // pragma omp parallel


    // compute inverse fft
    ierr=VecGetArray(xsmooth,&p_xsmooth); CHKERRQ(ierr);
    accfft_execute_c2r_t<ScalarTypeFD,ScalarType>(this->m_Opt->GetFFT().plan,this->m_xhat,p_xsmooth,timer);
    ierr=VecRestoreArray(xsmooth,&p_xsmooth); CHKERRQ(ierr);

    // increment fft timer
    this->m_Opt->IncreaseFFTTimers(timer);

    PetscFunctionReturn(0);

}




} // end of namespace

#endif // _PREPROCREG_CPP_
