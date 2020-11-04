#include <vector>
#include <string>

namespace reg {
  class RegOpt;
  class CLAIREInterface;
  class ReadWriteReg;
}

struct data_container;

class claire {
public:
  claire();
  ~claire();
  
  void setParameters(const std::vector<std::string>& args);
  
  void setReferenceImage(const double* data, size_t size);
  void setTemplateImage(const double* data, size_t size);
  void setVelocityField(const double* data, size_t size);
  
  void getVelocityField(double* data, size_t size);
  void getFinalState(double* data, size_t size);
  
  void runRegistration();
  void runForwardSolver();
private:
  reg::RegOpt* regopt;
  reg::CLAIREInterface* registration;
  reg::ReadWriteReg* readwrite;
  data_container* fields;
};
