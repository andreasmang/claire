#include "PythonInterface.hpp"
#include "RegOpt.hpp"
#include "CLAIREUtils.hpp"
#include "CLAIREInterface.hpp"
#include "ReadWriteReg.hpp"

struct data_container {
  data_container() { mR = nullptr; mT = nullptr; v = nullptr; aux = nullptr; v1 = nullptr; v2 = nullptr; v3 = nullptr; }
  ~data_container() {
    if (mR) VecDestroy(&mR);
    if (mT) VecDestroy(&mT);
    if (v1) VecDestroy(&v1);
    if (v2) VecDestroy(&v2);
    if (v3) VecDestroy(&v3);
    if (v) delete v;
    if (aux) VecDestroy(&aux);
  }
  Vec mR;
  Vec mT;
  Vec v1;
  Vec v2;
  Vec v3;
  reg::VecField* v;
  Vec aux;
};

void CopyFrom(double* dst, Vec src, size_t size) {
  const ScalarType* ptr;
  
  VecGetArrayRead(src, &ptr);

  if (sizeof(double) == sizeof(ScalarType)) {
    std::copy(ptr, ptr+size, dst);
  } else {
    for (size_t i=0;i<size;++i) dst[i] = static_cast<double>(ptr[i]);
  }

  VecRestoreArrayRead(src, &ptr);
}

void CopyTo(Vec dst, const double* src, size_t size) {
  ScalarType* ptr;
  
  VecGetArray(dst, &ptr);

  if (sizeof(double) == sizeof(ScalarType)) {
    std::copy(src, src+size, ptr);
  } else {
    for (size_t i=0;i<size;++i) ptr[i] = static_cast<double>(src[i]);
  }

  VecRestoreArray(dst, &ptr);
}

static size_t claire_cnt = 0;

claire::claire() {
  if (claire_cnt == 0) {
    PetscInitialize(0, reinterpret_cast<char***>(NULL),
                       reinterpret_cast<char*>(NULL),
                       reinterpret_cast<char*>(NULL));
  }
  claire_cnt++;
  
  regopt = new reg::RegOpt();
  readwrite = new reg::ReadWriteReg(regopt);
  registration = new reg::CLAIREInterface(regopt);
  registration->SetReadWrite(readwrite);
  
  fields = new data_container();  
}
claire::~claire() {
  delete regopt;
  delete readwrite;
  delete registration;
  delete fields;
  
  claire_cnt--;
  if (claire_cnt == 0) {
    reg::Finalize();
  }
}

void claire::setParameters(const std::vector<std::string>& args) {
  regopt->SetParameter(args);
  if (!regopt->m_SetupDone) {
    regopt->DoSetup();
  }
  if (!fields->v) {
    reg::VecCreate(fields->v1, regopt->m_Domain.nl, regopt->m_Domain.ng);
    reg::VecCreate(fields->v2, regopt->m_Domain.nl, regopt->m_Domain.ng);
    reg::VecCreate(fields->v3, regopt->m_Domain.nl, regopt->m_Domain.ng);
    fields->v = new reg::VecField(regopt, fields->v1, fields->v2, fields->v3);
    fields->v->SetValue(0.);
  }
  registration->SetInitialGuess(fields->v, false);
  if (!fields->mT) reg::VecCreate(fields->mT, regopt->m_Domain.nl, regopt->m_Domain.ng);
  registration->SetTemplateImage(fields->mT);
  if (!fields->mR) reg::VecCreate(fields->mR, regopt->m_Domain.nl, regopt->m_Domain.ng);
  registration->SetReferenceImage(fields->mR);
}

void claire::setReferenceImage(const double* data, size_t size) {
  if (size != static_cast<size_t>(regopt->m_Domain.nl)) throw;
  CopyTo(fields->mR, data, size);
  registration->SetReferenceImage(fields->mR);
}
void claire::setTemplateImage(const double* data, size_t size) {
  if (size != static_cast<size_t>(regopt->m_Domain.nl)) throw;
  CopyTo(fields->mT, data, size);
  registration->SetTemplateImage(fields->mT);
}
void claire::setVelocityField(const double* data, size_t size) {
  if (size != static_cast<size_t>(regopt->m_Domain.nl)*3) throw;
  CopyTo(fields->v1, data, regopt->m_Domain.nl);
  CopyTo(fields->v2, data+regopt->m_Domain.nl, regopt->m_Domain.nl);
  CopyTo(fields->v3, data+2*regopt->m_Domain.nl, regopt->m_Domain.nl);
  registration->SetInitialGuess(fields->v, false);
}

void claire::getVelocityField(double* data, size_t size) {
  if (size != static_cast<size_t>(regopt->m_Domain.nl)*3) throw;
  CopyFrom(data, fields->v1, regopt->m_Domain.nl);
  CopyFrom(data+regopt->m_Domain.nl, fields->v2, regopt->m_Domain.nl);
  CopyFrom(data+2*regopt->m_Domain.nl, fields->v3, regopt->m_Domain.nl);
}
void claire::getFinalState(double* data, size_t size) {
  if (size != static_cast<size_t>(regopt->m_Domain.nl)) throw;
  if (!fields->aux) reg::VecCreate(fields->aux, regopt->m_Domain.nl, regopt->m_Domain.ng);
  registration->GetFinalState(fields->aux);
  CopyFrom(data, fields->aux, size);
}

void claire::runRegistration() {
  registration->Run();
}
void claire::runForwardSolver() {
  registration->SolveForwardProblem();
}
