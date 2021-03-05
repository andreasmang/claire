#include "TwoLevel.hpp"

template<typename T>
__global__ void Restrict0stCentralGPU(T* __restrict__ dst,
                                      const T* __restrict__ src,
                                      dim3 nx) {
  T val;
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_x[2], idx_y[2], idx_z[2];
  
  idx_x[0] = gidx.x*2;
  idx_x[1] = idx_x[0] + 1;
  
  idx_y[0] = gidx.y*2;
  idx_y[1] = idx_y[0] + 1;
  
  idx_z[0] = gidx.z*2;
  idx_z[1] = idx_z[0] + 1;
  
  val = 0;
  
  for (int x=0; x<2; ++x) {
    unsigned int base_x = idx_x[x]*stride_x;
    for (int y=0; y<2; ++y) {
      const T* ptr = &src[base_x + idx_y[y]*stride_y];
      val += ptr[idx_z[0]] + ptr[idx_z[1]];
    }
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = 0.125*val;
}

template<typename T>
__global__ void Restrict1stIncludeGPU(T* __restrict__ dst,
                                        const T* __restrict__ src,
                                        dim3 nx) {
  T stage0[3];
  T stage1[3];
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_x[3], idx_y[3], idx_z[3];
  
  idx_x[0] = (gidx.x == 0 ? nx.x - 1 : gidx.x*2 - 1);
  idx_x[1] = gidx.x*2;
  idx_x[2] = idx_x[1] + 1;
  
  idx_y[0] = (gidx.y == 0 ? nx.y - 1 : gidx.y*2 - 1);
  idx_y[1] = gidx.y*2;
  idx_y[2] = idx_y[1] + 1;
  
  idx_z[0] = (gidx.z == 0 ? nx.z - 1 : gidx.z*2 - 1);
  idx_z[1] = gidx.z*2;
  idx_z[2] = idx_z[1] + 1;
  
  for (int x=0; x<3; ++x) {
    unsigned int base_x = idx_x[x]*stride_x;
    for (int y=0; y<3; ++y) {
      const T* ptr = &src[base_x + idx_y[y]*stride_y];
      stage0[y] = 0.25*(ptr[idx_z[0]] + 2.*ptr[idx_z[1]] + ptr[idx_z[2]]);
    }
    stage1[x] = 0.25*(stage0[0] + 2.*stage0[1] + stage0[2]);
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = 0.25*(stage1[0] + 2.*stage1[1] + stage1[2]);
}

template<typename T>
__global__ void Restrict1stCentralGPU(T* __restrict__ dst,
                                        const T* __restrict__ src,
                                        dim3 nx) {
  T stage0[4];
  T stage1[4];
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_x[4], idx_y[4], idx_z[4];
  
  idx_x[0] = (gidx.x == 0 ? nx.x - 1 : gidx.x*2 - 1);
  idx_x[1] = gidx.x*2;
  idx_x[2] = idx_x[1] + 1;
  idx_x[3] = (idx_x[2] + 1 == nx.x ? 0 : idx_x[2] + 1);
  
  idx_y[0] = (gidx.y == 0 ? nx.y - 1 : gidx.y*2 - 1);
  idx_y[1] = gidx.y*2;
  idx_y[2] = idx_y[1] + 1;
  idx_y[3] = (idx_y[2] + 1 == nx.y ? 0 : idx_y[2] + 1);
  
  idx_z[0] = (gidx.z == 0 ? nx.z - 1 : gidx.z*2 - 1);
  idx_z[1] = gidx.z*2;
  idx_z[2] = idx_z[1] + 1;
  idx_z[3] = (idx_z[2] + 1 == nx.z ? 0 : idx_z[2] + 1);
  
  for (int x=0; x<4; ++x) {
    unsigned int base_x = idx_x[x]*stride_x;
    for (int y=0; y<4; ++y) {
      const T* ptr = &src[base_x + idx_y[y]*stride_y];
      stage0[y] = 0.125*(ptr[idx_z[0]] + 3.*ptr[idx_z[1]] + 3.*ptr[idx_z[2]] + ptr[idx_z[3]]);
    }
    stage1[x] = 0.125*(stage0[0] + 3.*stage0[1] + 3.*stage0[2] + stage0[3]);
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = 0.125*(stage1[0] + 3.*stage1[1] + 3.*stage1[2] + stage1[3]);
}

template<typename T>
__global__ void Restrict2ndCentralGPU(T* __restrict__ dst,
                                      const T* __restrict__ src,
                                      dim3 nx) {
  T stage0[6];
  T stage1[6];
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_x[6], idx_y[6], idx_z[6];
  
  idx_x[0] = (gidx.x == 0 ? nx.x - 2 : gidx.x*2 - 2);
  idx_x[1] = (gidx.x == 0 ? nx.x - 1 : gidx.x*2 - 1);
  idx_x[2] = gidx.x*2;
  idx_x[3] = idx_x[2] + 1;
  idx_x[4] = (idx_x[3] + 1 == nx.x ? 0 : idx_x[3] + 1);
  idx_x[5] = (idx_x[4] + 1 == nx.x ? 0 : idx_x[4] + 1);
  
  idx_y[0] = (gidx.y == 0 ? nx.y - 2 : gidx.y*2 - 2);
  idx_y[1] = (gidx.y == 0 ? nx.y - 1 : gidx.y*2 - 1);
  idx_y[2] = gidx.y*2;
  idx_y[3] = idx_y[2] + 1;
  idx_y[4] = (idx_y[3] + 1 == nx.y ? 0 : idx_y[3] + 1);
  idx_y[5] = (idx_y[4] + 1 == nx.y ? 0 : idx_y[4] + 1);
  
  idx_z[0] = (gidx.z == 0 ? nx.z - 2 : gidx.z*2 - 2);
  idx_z[1] = (gidx.z == 0 ? nx.z - 1 : gidx.z*2 - 1);
  idx_z[2] = gidx.z*2;
  idx_z[3] = idx_z[2] + 1;
  idx_z[4] = (idx_z[3] + 1 == nx.z ? 0 : idx_z[3] + 1);
  idx_z[5] = (idx_z[4] + 1 == nx.z ? 0 : idx_z[4] + 1);
  
  for (int x=0; x<6; ++x) {
    unsigned int base_x = idx_x[x]*stride_x;
    for (int y=0; y<6; ++y) {
      const T* ptr = &src[base_x + idx_y[y]*stride_y];
      stage0[y] = (-3.*ptr[idx_z[0]]+5.*ptr[idx_z[1]]+30.*ptr[idx_z[2]]+30.*ptr[idx_z[3]]+5.*ptr[idx_z[4]]-3.*ptr[idx_z[5]])/64.;
    }
    stage1[x] = (-3.*stage0[0]+5.*stage0[1]+30.*stage0[2]+30.*stage0[3]+5.*stage0[4]-3.*stage0[5])/64.;
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = (-3.*stage1[0]+5.*stage1[1]+30.*stage1[2]+30.*stage1[3]+5.*stage1[4]-3.*stage1[5])/64.;
}

template<typename T>
__global__ void Restrict3rdCentralGPU(T* __restrict__ dst,
                                      const T* __restrict__ src,
                                      dim3 nx) {
  T stage0[8];
  T stage1[8];
  
  dim3 gidx;
  
  gidx.x = blockIdx.z;
  gidx.y = blockIdx.y;
  gidx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (gidx.z*2 >= nx.z) return;
  
  unsigned int stride_y = nx.z;
  unsigned int stride_x = nx.z*nx.y;
  
  unsigned int idx_z[8];
  
  idx_z[0] = (gidx.z*2 + nx.z - 3)%nx.z;
  idx_z[1] = (gidx.z*2 + nx.z - 2)%nx.z;
  idx_z[2] = (gidx.z*2 + nx.z - 1)%nx.z;
  idx_z[3] = (gidx.z*2);
  idx_z[4] = (gidx.z*2 + 1);
  idx_z[5] = (gidx.z*2 + 2)%nx.z;
  idx_z[6] = (gidx.z*2 + 3)%nx.z;
  idx_z[7] = (gidx.z*2 + 4)%nx.z;
  
  for (int x=0; x<8; ++x) {
    unsigned int idx_x = (gidx.x*2 + x + nx.x - 3)%nx.x;
    unsigned int base_x = idx_x*stride_x;
    for (int y=0; y<8; ++y) {
      unsigned int idx_y = (gidx.y*2 + y + nx.y - 3)%nx.y;
      const T* ptr = &src[base_x + idx_y*stride_y];
      stage0[y] = (
        -3.  *ptr[idx_z[0]]
        -9.  *ptr[idx_z[1]]
        +29. *ptr[idx_z[2]]
        +111.*ptr[idx_z[3]]
        +111.*ptr[idx_z[4]]
        +29. *ptr[idx_z[5]]
        -9.  *ptr[idx_z[6]]
        -3.  *ptr[idx_z[7]]
      )/256.;
    }
    stage1[x] = (
        -3.  *stage0[0]
        -9.  *stage0[1]
        +29. *stage0[2]
        +111.*stage0[3]
        +111.*stage0[4]
        +29. *stage0[5]
        -9.  *stage0[6]
        -3.  *stage0[7]
      )/256.;
  }
  
  stride_y /= 2;
  stride_x /= 4;
  dst[gidx.z + gidx.y*stride_y + gidx.x*stride_x] = (
        -3.  *stage1[0]
        -9.  *stage1[1]
        +29. *stage1[2]
        +111.*stage1[3]
        +111.*stage1[4]
        +29. *stage1[5]
        -9.  *stage1[6]
        -3.  *stage1[7]
      )/256.;
}

template<typename T>
__global__ void Prolong0stCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val;
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  val = src[idx.x*stride_x + idx.y*stride_y + idx.z];

  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val;
  dst[base + 1 + 0*stride_y + 0*stride_x] = val;
  dst[base + 0 + 1*stride_y + 0*stride_x] = val;
  dst[base + 1 + 1*stride_y + 0*stride_x] = val;
  dst[base + 0 + 0*stride_y + 1*stride_x] = val;
  dst[base + 1 + 0*stride_y + 1*stride_x] = val;
  dst[base + 0 + 1*stride_y + 1*stride_x] = val;
  dst[base + 1 + 1*stride_y + 1*stride_x] = val;
}

template<typename T>
__global__ void Prolong1stIncludeGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  unsigned int idx_z[2], idx_y[2], idx_x[2];
  
  idx_z[0] =  idx.z;
  idx_z[1] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);
  
  idx_y[0] =  idx.y*stride_y;
  idx_y[1] = (idx.y + 1 < nx.y ? idx.y + 1 : 0)*stride_y;
  
  idx_x[0] =  idx.x*stride_x;
  idx_x[1] = (idx.x + 1 < nx.x ? idx.x + 1 : 0)*stride_x;

  val[0] = src[idx_z[0] + idx_y[0] + idx_x[0]];
  val[1] = src[idx_z[1] + idx_y[0] + idx_x[0]];
  val[2] = src[idx_z[0] + idx_y[1] + idx_x[0]];
  val[3] = src[idx_z[1] + idx_y[1] + idx_x[0]];
  val[4] = src[idx_z[0] + idx_y[0] + idx_x[1]];
  val[5] = src[idx_z[1] + idx_y[0] + idx_x[1]];
  val[6] = src[idx_z[0] + idx_y[1] + idx_x[1]];
  val[7] = src[idx_z[1] + idx_y[1] + idx_x[1]];
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = (val[1] + val[0])*0.5;
  dst[base + 0 + 1*stride_y + 0*stride_x] = (val[2] + val[0])*0.5;
  dst[base + 1 + 1*stride_y + 0*stride_x] = (val[3] + val[2] + val[1] + val[0])*0.25;
  dst[base + 0 + 0*stride_y + 1*stride_x] = (val[4] + val[0])*0.5;
  dst[base + 1 + 0*stride_y + 1*stride_x] = (val[5] + val[4] + val[1] + val[0])*0.25;
  dst[base + 0 + 1*stride_y + 1*stride_x] = (val[6] + val[4] + val[2] + val[0])*0.25;
  dst[base + 1 + 1*stride_y + 1*stride_x] = (val[7] + val[6] + val[5] + val[4] + val[3] + val[2] + val[1] + val[0])*0.125;
}

template<typename T>
__global__ void Prolong1stCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  const T weight[3] = {0.25, 0.75, 0.};
  val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.;
  
  unsigned int idx_z[3];
  
  idx_z[0] = (idx.z >= 1 ? idx.z - 1 : nx.z - 1);
  idx_z[1] =  idx.z;
  idx_z[2] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);

  for (int x = 0; x < 3; ++x) {
    unsigned int idx_x = idx.x + x;
    idx_x = (idx_x >= 1 ? (idx_x - 1 >= nx.x ? idx_x - 1 - nx.x : idx_x - 1) : nx.x - 1 + idx_x);
    idx_x *= stride_x;
    for (int y = 0; y < 3; ++y) {
      unsigned int idx_y = idx.y + y;
      idx_y = (idx_y >= 1 ? (idx_y - 1 >= nx.y ? idx_y - 1 - nx.y : idx_y - 1) : nx.y - 1 + idx_y);
      const T* ptr = &src[idx_x + idx_y*stride_y];
      T cw[4];
      cw[0] = weight[y]  *weight[x];
      cw[1] = weight[2-y]*weight[x];
      cw[2] = weight[y]  *weight[2-x];
      cw[3] = weight[2-y]*weight[2-x];
      for (int z = 0; z < 3; ++z) {
        T inval = ptr[idx_z[z]];
        val[0] += inval*weight[z]  *cw[0];
        val[1] += inval*weight[2-z]*cw[0];
        val[2] += inval*weight[z]  *cw[1];
        val[3] += inval*weight[2-z]*cw[1];
        val[4] += inval*weight[z]  *cw[2];
        val[5] += inval*weight[2-z]*cw[2];
        val[6] += inval*weight[z]  *cw[3];
        val[7] += inval*weight[2-z]*cw[3];
      }
    }
  }
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = val[1];
  dst[base + 0 + 1*stride_y + 0*stride_x] = val[2];
  dst[base + 1 + 1*stride_y + 0*stride_x] = val[3];
  dst[base + 0 + 0*stride_y + 1*stride_x] = val[4];
  dst[base + 1 + 0*stride_y + 1*stride_x] = val[5];
  dst[base + 0 + 1*stride_y + 1*stride_x] = val[6];
  dst[base + 1 + 1*stride_y + 1*stride_x] = val[7];
}

template<typename T>
__global__ void Prolong2ndCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  const T weight[3] = {5./32., 30./32., -3./32.};
  val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.;
  
  unsigned int idx_z[3];
  
  idx_z[0] = (idx.z >= 1 ? idx.z - 1 : nx.z - 1);
  idx_z[1] =  idx.z;
  idx_z[2] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);

  for (int x = 0; x < 3; ++x) {
    unsigned int idx_x = idx.x + x;
    idx_x = (idx_x >= 1 ? (idx_x - 1 >= nx.x ? idx_x - 1 - nx.x : idx_x - 1) : nx.x - 1 + idx_x);
    idx_x *= stride_x;
    for (int y = 0; y < 3; ++y) {
      unsigned int idx_y = idx.y + y;
      idx_y = (idx_y >= 1 ? (idx_y - 1 >= nx.y ? idx_y - 1 - nx.y : idx_y - 1) : nx.y - 1 + idx_y);
      const T* ptr = &src[idx_x + idx_y*stride_y];
      T cw[4];
      cw[0] = weight[y]  *weight[x];
      cw[1] = weight[2-y]*weight[x];
      cw[2] = weight[y]  *weight[2-x];
      cw[3] = weight[2-y]*weight[2-x];
      for (int z = 0; z < 3; ++z) {
        T inval = ptr[idx_z[z]];
        val[0] += inval*weight[z]  *cw[0];
        val[1] += inval*weight[2-z]*cw[0];
        val[2] += inval*weight[z]  *cw[1];
        val[3] += inval*weight[2-z]*cw[1];
        val[4] += inval*weight[z]  *cw[2];
        val[5] += inval*weight[2-z]*cw[2];
        val[6] += inval*weight[z]  *cw[3];
        val[7] += inval*weight[2-z]*cw[3];
      }
    }
  }
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = val[1];
  dst[base + 0 + 1*stride_y + 0*stride_x] = val[2];
  dst[base + 1 + 1*stride_y + 0*stride_x] = val[3];
  dst[base + 0 + 0*stride_y + 1*stride_x] = val[4];
  dst[base + 1 + 0*stride_y + 1*stride_x] = val[5];
  dst[base + 0 + 1*stride_y + 1*stride_x] = val[6];
  dst[base + 1 + 1*stride_y + 1*stride_x] = val[7];
}


template<typename T>
__global__ void Prolong3rdCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  const T weight[5] = {-3./128., 29./128., 111./128., -9./128., 0.};
  val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.;
  
  unsigned int idx_z[5];
  
  idx_z[0] = (idx.z >= 2 ? idx.z - 2 : nx.z - 2 + idx.z);
  idx_z[1] = (idx.z >= 1 ? idx.z - 1 : nx.z - 1);
  idx_z[2] =  idx.z;
  idx_z[3] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);
  idx_z[4] = (idx.z + 2 < nx.z ? idx.z + 2 : idx.z + 2 - nx.z);

  for (int x = 0; x < 5; ++x) {
    unsigned int idx_x = idx.x + x;
    idx_x = (idx_x >= 2 ? (idx_x - 2 >= nx.x ? idx_x - 2 - nx.x : idx_x - 2) : nx.x - 2 + idx_x);
    idx_x *= stride_x;
    for (int y = 0; y < 5; ++y) {
      unsigned int idx_y = idx.y + y;
      idx_y = (idx_y >= 2 ? (idx_y - 2 >= nx.y ? idx_y - 2 - nx.y : idx_y - 2) : nx.y - 2 + idx_y);
      const T* ptr = &src[idx_x + idx_y*stride_y];
      T cw[4];
      cw[0] = weight[y]  *weight[x];
      cw[1] = weight[4-y]*weight[x];
      cw[2] = weight[y]  *weight[4-x];
      cw[3] = weight[4-y]*weight[4-x];
      for (int z = 0; z < 5; ++z) {
        T inval = ptr[idx_z[z]];
        val[0] += inval*weight[z]  *cw[0];
        val[1] += inval*weight[4-z]*cw[0];
        val[2] += inval*weight[z]  *cw[1];
        val[3] += inval*weight[4-z]*cw[1];
        val[4] += inval*weight[z]  *cw[2];
        val[5] += inval*weight[4-z]*cw[2];
        val[6] += inval*weight[z]  *cw[3];
        val[7] += inval*weight[4-z]*cw[3];
      }
    }
  }
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = val[1];
  dst[base + 0 + 1*stride_y + 0*stride_x] = val[2];
  dst[base + 1 + 1*stride_y + 0*stride_x] = val[3];
  dst[base + 0 + 0*stride_y + 1*stride_x] = val[4];
  dst[base + 1 + 0*stride_y + 1*stride_x] = val[5];
  dst[base + 0 + 1*stride_y + 1*stride_x] = val[6];
  dst[base + 1 + 1*stride_y + 1*stride_x] = val[7];
}

template<typename T>
__global__ void Prolong5thCentralGPU(T* __restrict__ dst,
                                       const T* __restrict__ src,
                                       dim3 nx) {
  T val[8];
  
  dim3 idx;
  
  idx.x = blockIdx.z;
  idx.y = blockIdx.y;
  idx.z = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (idx.z >= nx.z) return;
  
  unsigned int stride_x = nx.y*nx.z;
  unsigned int stride_y = nx.z;
  
  const T weight[7] = {63./8192., -495./8192., 2310./8192., 6930./8192., -693./8192., 77./8192., 0.};
  val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.;
  
  unsigned int idx_z[7];
  
  idx_z[0] = (idx.z >= 3 ? idx.z - 3 : nx.z - 3 + idx.z);
  idx_z[1] = (idx.z >= 2 ? idx.z - 2 : nx.z - 2 + idx.z);
  idx_z[2] = (idx.z >= 1 ? idx.z - 1 : nx.z - 1);
  idx_z[3] =  idx.z;
  idx_z[4] = (idx.z + 1 < nx.z ? idx.z + 1 : 0);
  idx_z[5] = (idx.z + 2 < nx.z ? idx.z + 2 : idx.z + 2 - nx.z);
  idx_z[6] = (idx.z + 3 < nx.z ? idx.z + 3 : idx.z + 3 - nx.z);

  for (int x = 0; x < 7; ++x) {
    unsigned int idx_x = idx.x + x;
    idx_x = (idx_x >= 3 ? (idx_x - 3 >= nx.x ? idx_x - 3 - nx.x : idx_x - 3) : nx.x - 3 + idx_x);
    idx_x *= stride_x;
    for (int y = 0; y < 7; ++y) {
      unsigned int idx_y = idx.y + y;
      idx_y = (idx_y >= 3 ? (idx_y - 3 >= nx.y ? idx_y - 3 - nx.y : idx_y - 3) : nx.y - 3 + idx_y);
      const T* ptr = &src[idx_x + idx_y*stride_y];
      T cw[4];
      cw[0] = weight[y]  *weight[x];
      cw[1] = weight[6-y]*weight[x];
      cw[2] = weight[y]  *weight[6-x];
      cw[3] = weight[6-y]*weight[6-x];
      for (int z = 0; z < 7; ++z) {
        T inval = ptr[idx_z[z]];
        val[0] += inval*weight[z]  *cw[0];
        val[1] += inval*weight[6-z]*cw[0];
        val[2] += inval*weight[z]  *cw[1];
        val[3] += inval*weight[6-z]*cw[1];
        val[4] += inval*weight[z]  *cw[2];
        val[5] += inval*weight[6-z]*cw[2];
        val[6] += inval*weight[z]  *cw[3];
        val[7] += inval*weight[6-z]*cw[3];
      }
    }
  }
  
  stride_y *= 2;
  stride_x *= 4;
  
  unsigned int base = (idx.z + idx.y*stride_y + idx.x*stride_x)*2;
  
  dst[base + 0 + 0*stride_y + 0*stride_x] = val[0];
  dst[base + 1 + 0*stride_y + 0*stride_x] = val[1];
  dst[base + 0 + 1*stride_y + 0*stride_x] = val[2];
  dst[base + 1 + 1*stride_y + 0*stride_x] = val[3];
  dst[base + 0 + 0*stride_y + 1*stride_x] = val[4];
  dst[base + 1 + 0*stride_y + 1*stride_x] = val[5];
  dst[base + 0 + 1*stride_y + 1*stride_x] = val[6];
  dst[base + 1 + 1*stride_y + 1*stride_x] = val[7];
}

namespace reg {
  
namespace TwoLevelFiniteKernel {

PetscErrorCode InitConstants(IntType* isize, IntType* isize_g, ScalarType* hx, IntType* halo) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode Restrict(ScalarType* dst, const ScalarType *src, IntType *nl_f) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  IntType nl[3];
  nl[0] = nl_f[0]/2; nl[1] = nl_f[1]/2; nl[2] = nl_f[2]/2;
  
  dim3 block, grid;
  
  if (nl[2] <= 1024 && nl[2] >= 32) {
    block.x = (nl[2] + 31)/32;
    block.x *= 32;
    grid.x = 1;
  } else {
    block.x = 128; // 128 threads per block
    grid.x = (nl[2] + 127)/128;  // $\lceil nl_2 / 128 \rceil$
  }
  grid.y = nl[1];
  grid.z = nl[0];
  dim3 nl3;
  nl3.x = nl_f[0]; nl3.y = nl_f[1]; nl3.z = nl_f[2];

  Restrict3rdCentralGPU<<<grid, block>>>(dst, src, nl3);
  ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
  ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);  
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode Prolong(ScalarType* dst, const ScalarType *src, IntType *nl_f) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  
  IntType nl[3];
  nl[0] = nl_f[0]/2; nl[1] = nl_f[1]/2; nl[2] = nl_f[2]/2;
  
  dim3 block, grid;
  
  if (nl[2] <= 1024 && nl[2] >= 32) {
    block.x = (nl[2] + 31)/32;
    block.x *= 32;
    grid.x = 1;
  } else {
    block.x = 128; // 128 threads per block
    grid.x = (nl[2] + 127)/128;  // $\lceil nl_2 / 128 \rceil$
  }
  grid.y = nl[1];
  grid.z = nl[0];
  dim3 nl3;
  nl3.x = nl[0]; nl3.y = nl[1]; nl3.z = nl[2];

  Prolong3rdCentralGPU<<<grid, block>>>(dst, src, nl3);
  ierr = cudaDeviceSynchronize(); CHKERRCUDA(ierr);
  ierr = cudaCheckKernelError(); CHKERRCUDA(ierr);  
  
  PetscFunctionReturn(ierr);
}

} // namespace TwoLevelFiniteKernel
} // namespace reg
