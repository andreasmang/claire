#include <interp3_gpu_mpi.hpp>

template <typename T>
struct square
{
    __host__ __device__
        T operator()(const T& x) const { 
            return x * x;
        }
};

void get_count(int * arr, int size, int val, int* count) {
  thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast<int>(arr);
  *count = thrust::count(dev_ptr, dev_ptr + size, val);
}

void print_norm(ScalarType* arr, int N) {
  thrust::device_ptr<ScalarType> d_arr = thrust::device_pointer_cast<ScalarType>(arr);
  square<ScalarType>  unary_op;
  thrust::plus<ScalarType> binary_op;
  ScalarType init = 0; 
  ScalarType norm = std::sqrt( thrust::transform_reduce(d_arr, d_arr+N, unary_op, init, binary_op) );
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "norm = %f\n", norm);
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void print_max(ScalarType *arr, int N) {
  thrust::device_ptr<ScalarType> d_arr = thrust::device_pointer_cast<ScalarType>(arr);
  ScalarType result = *(thrust::max_element(thrust::device, d_arr, d_arr+N));
  PetscPrintf(PETSC_COMM_WORLD, "%f\n", result);
}


void test_count() {

  int dev_id;
  cudaGetDevice(&dev_id);
  int count = 0;
  thrust::device_ptr<int> void_ptr = thrust::device_malloc<int>(10);
//    std::cout << "rank " << dev_id << " allocated memory " << std::endl;
  
  thrust::fill(void_ptr, void_ptr + 10, 10);
  count = thrust::count(void_ptr, void_ptr + 10, 10);
  std::cout << "[" << dev_id << "] = " << count << std::endl;

  thrust::device_free(void_ptr);
//    std::cout << "rank " << dev_id << " freed  memory " << std::endl;

  //thrust::device_vector<int> vec(10, 2);
  //for (int i=0; i<10; i++) {
  //    std::cout <<  vec[i] << std::endl;
  //}
  //count = thrust::count(vec.begin(), vec.end(), 2);
}
