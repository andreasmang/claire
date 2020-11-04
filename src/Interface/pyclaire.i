%module pyclaire

%{
#define SWIG_FILE_WITH_INIT
#include "PythonInterface.hpp"
%}

%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(StringStringMap) vector<string>;
}

%include "numpy.i"
%init %{
import_array();
%}
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%apply (double* INPLACE_ARRAY_FLAT, size_t DIM_FLAT) {(const double *data, size_t size)};
%apply (double* INPLACE_ARRAY_FLAT, size_t DIM_FLAT) {(double *data, size_t size)};

%include "PythonInterface.hpp"
