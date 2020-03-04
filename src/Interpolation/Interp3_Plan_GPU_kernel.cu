#include <interp3_gpu_mpi.hpp>


void get_count(int * arr, int size, int val, int* count) {
    thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast<int>(arr);
    *count = thrust::count(dev_ptr, dev_ptr + size, val);
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
