// This function performs a 3D cubic interpolation.

#define _mm256_set_m128(va, vb) \
          _mm256_insertf128_ps(_mm256_castps128_ps256(vb), va, 1)
#define _mm512_set_m256(va, vb) \
          _mm512_insertf32x8(_mm512_castps256_ps512(vb), va, 1)

#include <cmath>
#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>

#include <interp3.hpp>
#ifdef FAST_INTERPV
#include <immintrin.h>
#endif
#define COORD_DIM 3
//#define VERBOSE2
#define sleep(x) ;

/*
 * multiply and shift query points based on each proc's domain
 */
void rescale_xyz(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		int* isize, int* isize_g, const int N_pts, Real* Q_) {

	if (g_size == 0)
		return;
	Real h[3];

	h[0] = 1. / (N_reg[0]);
	h[1] = 1. / (N_reg[1]);
	h[2] = 1. / (N_reg[2]);

  const Real iX0 = istart[0]*h[0];
  const Real iX1 = istart[1]*h[1];
  const Real iX2 = istart[2]*h[2];
  const Real N0 = N_reg[0];
  const Real N1 = N_reg[1];
  const Real N2 = N_reg[2];

#pragma omp parallel for
	for (int i = 0; i < N_pts; i++) {
    Real* Q_ptr = &Q_[COORD_DIM * i];
    Q_ptr[0] = (Q_ptr[0]-iX0)*N0+g_size;
    Q_ptr[1] = (Q_ptr[1]-iX1)*N1+g_size;
    Q_ptr[2] = (Q_ptr[2]-iX2)*N2+g_size;
	}

  //std::cout <<std::floor(N_pts/16.0)*16<< '\t' <<std::ceil(N_pts/16.0)*16 << std::endl;
  // set Q to be one so that for the peeled loop we only access grid indxx =0
  if(N_pts%16 != 0)
	for (int i = std::floor(N_pts/16.0)*16+N_pts%16; i < std::ceil(N_pts/16.0)*16; i++) {
    Real* Q_ptr = &Q_[COORD_DIM * i];
    Q_ptr[0] = 1;
    Q_ptr[1] = 1;
    Q_ptr[2] = 1;
	}
	return;
} // end of rescale_xyz

void rescale_xyzgrid(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		int* isize, int* isize_g, const int N_pts, pvfmm::Iterator<Real> Q_) {

	if (g_size == 0)
		return;
	Real h[3];

	h[0] = 1. / (N_reg[0]);
	h[1] = 1. / (N_reg[1]);
	h[2] = 1. / (N_reg[2]);

  const Real iX0 = istart[0]*h[0];
  const Real iX1 = istart[1]*h[1];
  const Real iX2 = istart[2]*h[2];
  const Real Nx = N_reg[0];
  const Real Ny = N_reg[1];
  const Real Nz = N_reg[2];
  const int isize_g2 = isize_g[2];
  const int isize_g1g2 = isize_g2 * isize_g[1];
  pvfmm::Iterator<Real> tmp = pvfmm::aligned_new<Real>(N_pts*3);
  pvfmm::memcopy(tmp, Q_, N_pts*3);

#pragma omp parallel for
	for (int i = 0; i < N_pts; i++) {
    Real* Q_ptr = &Q_[4 * i];
    Real* tmp_ptr = &tmp[3 * i];
    // std::cout << tmp_ptr[0] << '\t' << tmp_ptr[1] << '\t' << tmp_ptr[2] << std::endl;
    Q_ptr[0] = (tmp_ptr[0]-iX0)*Nx+g_size;
    Q_ptr[1] = (tmp_ptr[1]-iX1)*Ny+g_size;
    Q_ptr[2] = (tmp_ptr[2]-iX2)*Nz+g_size;

		const int grid_indx0 = ((int)(Q_ptr[0])) - 1;
		Q_ptr[0] -= grid_indx0;
		const int grid_indx1 = ((int)(Q_ptr[1])) - 1;
		Q_ptr[1] -= grid_indx1;
		const int grid_indx2 = ((int)(Q_ptr[2])) - 1;
		Q_ptr[2] -= grid_indx2;
		const int indxx = isize_g1g2 * grid_indx0 + grid_indx2 + isize_g2 * grid_indx1 ;
    Q_ptr[3] = (Real)indxx;
    // std::cout << grid_indx0 << '\t'
    // << grid_indx1 << '\t'
    // << grid_indx2 << std::endl;
	}
  pvfmm::aligned_delete(tmp);

  //std::cout <<std::floor(N_pts/16.0)*16<< '\t' <<std::ceil(N_pts/16.0)*16 << std::endl;
  // set Q to be one so that for the peeled loop we only access grid indxx =0
  if(N_pts%16 != 0)
	for (int i = std::floor(N_pts/16.0)*16+N_pts%16; i < std::ceil(N_pts/16.0)*16; i++) {
    Real* Q_ptr = &Q_[4 * i];
    Q_ptr[0] = 1;
    Q_ptr[1] = 1;
    Q_ptr[2] = 1;
    Q_ptr[3] = 0;
	}
	return;
} // end of rescale_xyz
// acknowledgemet to http://stackoverflow.com/questions/13219146/how-to-sum-m256-horizontally
// x = ( x7, x6, x5, x4, x3, x2, x1, x0 )
#ifdef FAST_INTERPV

float sum8(__m256 x) {
    // hiQuad = ( x7, x6, x5, x4 )
    const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
    // loQuad = ( x3, x2, x1, x0 )
    const __m128 loQuad = _mm256_castps256_ps128(x);
    // sumQuad = ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
    const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
    // loDual = ( -, -, x1 + x5, x0 + x4 )
    const __m128 loDual = sumQuad;
    // hiDual = ( -, -, x3 + x7, x2 + x6 )
    const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
    // sumDual = ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
    const __m128 sumDual = _mm_add_ps(loDual, hiDual);
    // lo = ( -, -, -, x0 + x2 + x4 + x6 )
    const __m128 lo = sumDual;
    // hi = ( -, -, -, x1 + x3 + x5 + x7 )
    const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
    // sum = ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
    const __m128 sum = _mm_add_ss(lo, hi);
    return _mm_cvtss_f32(sum);
}
void print128(__m128 x, const char* name) {
  Real* ptr = (Real*)&x;
  std::cout << name
            << " [0] = " << ptr[0]
            << " [1] = " << ptr[1]
            << " [2] = " << ptr[2]
            << " [3] = " << ptr[3] << std::endl;
}
void print256(__m256 x, const char* name) {
  Real* ptr = (Real*)&x;
  std::cout << name
            << "\n [0] = " << ptr[0]
            << "\n [1] = " << ptr[1]
            << "\n [2] = " << ptr[2]
            << "\n [3] = " << ptr[3]
            << "\n [4] = " << ptr[4]
            << "\n [5] = " << ptr[5]
            << "\n [6] = " << ptr[6]
            << "\n [7] = " << ptr[7] << std::endl;
}
void print512(__m512 &x, const char* name) {
  Real* ptr = (Real*)&x;
  std::cout << name
            << "\n [0] = " << ptr[0]
            << "\n [1] = " << ptr[1]
            << "\n [2] = " << ptr[2]
            << "\n [3] = " << ptr[3]
            << "\n [4] = " << ptr[4]
            << "\n [5] = " << ptr[5]
            << "\n [6] = " << ptr[6]
            << "\n [7] = " << ptr[7]
            << "\n [8] = " << ptr[8]
            << "\n [9] = " << ptr[9]
            << "\n [10] = " << ptr[10]
            << "\n [11] = " << ptr[11]
            << "\n [12] = " << ptr[12]
            << "\n [13] = " << ptr[13]
            << "\n [14] = " << ptr[14]
            << "\n [15] = " << ptr[15] << std::endl;
}

#endif



#ifdef FAST_INTERPV
#if defined(KNL)
void vectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled) {

#ifdef INTERP_DEBUG
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  PCOUT << "In KNL kernel\n";
#endif
  const __m512  c1000_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-1.0,-0.0,-0.0,-0.0));
  const __m512  c2211_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-2.0,-2.0,-1.0,-1.0));
  const __m512  c3332_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-3.0,-3.0,-3.0,-2.0));

  const __m512  c1000_512_ = _mm512_set_ps(
      -1.0,-1.0,-1.0,-1.0,
      -0.0,-0.0,-0.0,-0.0,
      -0.0,-0.0,-0.0,-0.0,
      -0.0,-0.0,-0.0,-0.0
      );
  const __m512  c2211_512_ = _mm512_set_ps(
      -2.0,-2.0,-2.0,-2.0,
      -2.0,-2.0,-2.0,-2.0,
      -1.0,-1.0,-1.0,-1.0,
      -1.0,-1.0,-1.0,-1.0
      );
  const __m512  c3332_512_ = _mm512_set_ps(
      -3.0,-3.0,-3.0,-3.0,
      -3.0,-3.0,-3.0,-3.0,
      -3.0,-3.0,-3.0,-3.0,
      -2.0,-2.0,-2.0,-2.0
      );
  const __m512 vlagr_0000_1111_2222_3333_512  = _mm512_set_ps(
      -0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,
       0.5,0.5,0.5,0.5,
      -0.5,-0.5,-0.5,-0.5,
      +0.1666666667,0.1666666667,0.1666666667,0.1666666667
      );
  //print512(c1000_512,"512");
  //print512(c3332_512,"");
  //do{}while(1);

  const __m512 vlagr_512 = _mm512_set_ps(
      -0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667,
      -0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667
      );
  const int isize_g2 = isize_g[2];
  const int two_isize_g2 = 2*isize_g2;
  const int three_isize_g2 = 3*isize_g2;
  const int reg_plus = isize_g[1]*isize_g2;
  const int NzNy = isize_g2 * isize_g[1];
  Real* Q_ptr = query_points;


  // std::cout << "KNL" << std::endl;
  //_mm_prefetch( (char*)Q_ptr,_MM_HINT_NTA);



  int CHUNK=16;

//#pragma omp parallel for
//	for (int ii = 0; ii < (int)std::ceil(N_pts/(float)CHUNK); ii++) {
#pragma omp parallel for
	for (int i = 0; i < N_pts; i++) {
#ifdef INTERP_USE_MORE_MEM_L1
		Real point[COORD_DIM];
		point[0] = Q_ptr[i*4+0];
		point[1] = Q_ptr[i*4+1];
		point[2] = Q_ptr[i*4+2];
    const int indxx = (int) Q_ptr[4*i + 3];
#else
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = Q_ptr[i*3+0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = Q_ptr[i*3+1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = Q_ptr[i*3+2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];
    // Q_ptr += 3;
		const int indxx = NzNy * grid_indx[0] + grid_indx[2] + isize_g2 * grid_indx[1] ;
#endif
    //_mm_prefetch( (char*)Q_ptr,_MM_HINT_T2);

    __m512 vM1_0000_1111_2222_3333(vlagr_0000_1111_2222_3333_512);
    __m512 vM2_512(vlagr_512);
    //__m512 vM2_512(vlagr_3333_2222_0000_1111_512);
    __m512 vM0_512(vlagr_512);

    {
    const __m512 vx0_512 =  _mm512_set1_ps(point[0]);
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c1000_512));
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c2211_512));
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c3332_512));

    const __m512 vx1_512 =  _mm512_set1_ps(point[1]);
    vM1_0000_1111_2222_3333  = _mm512_mul_ps(vM1_0000_1111_2222_3333 , _mm512_add_ps(vx1_512,c1000_512_));
    vM1_0000_1111_2222_3333  = _mm512_mul_ps(vM1_0000_1111_2222_3333 , _mm512_add_ps(vx1_512,c2211_512_));
    vM1_0000_1111_2222_3333  = _mm512_mul_ps(vM1_0000_1111_2222_3333 , _mm512_add_ps(vx1_512,c3332_512_));

    //const __m512 vx2_512 =  _mm512_set1_ps(point[2]);
    //vM2_512  = _mm512_mul_ps(vM2_512, _mm512_add_ps(vx2_512,c0010_512));
    //vM2_512  = _mm512_mul_ps(vM2_512, _mm512_add_ps(vx2_512,c1122_512));
    //vM2_512  = _mm512_mul_ps(vM2_512, _mm512_add_ps(vx2_512,c3233_512));
    //vM0_512  = _mm512_mul_ps(vM0_512 , vM2_512);
    const __m512 vx2_512 =  _mm512_set1_ps(point[2]);
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c1000_512));
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c2211_512));
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c3332_512));
    //vM0_512  = _mm512_mul_ps(vM0_512 , vM2_512);
    //print512(vM0_512, "vM0");
    //print512(vM2_512, "vM2");
    //print512(vM2_512, "vM2");
    //print512(vM1_0000_1111_2222_3333  , "vM1");
    //do{}while(1);
    }

    int indx = 0;
    Real* reg_ptr = reg_grid_vals + indxx;//&reg_grid_vals[indxx];
    //_mm_prefetch( (char*)reg_ptr,_MM_HINT_T0);

    // load all vfij

          __m512 vf_i0_j0123 = _mm512_setzero_ps();
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b1111000000000000, reg_ptr);
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000000011110000, reg_ptr+two_isize_g2);
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000000000001111, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          __m512 vf_i1_j0123 = _mm512_setzero_ps();
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b1111000000000000, reg_ptr);
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000000011110000, reg_ptr+two_isize_g2);
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000000000001111, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          __m512 vf_i2_j0123 = _mm512_setzero_ps();
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b1111000000000000, reg_ptr);
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000000011110000, reg_ptr+two_isize_g2);
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000000000001111, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          __m512 vf_i3_j0123 = _mm512_setzero_ps();
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b1111000000000000, reg_ptr);
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000000011110000, reg_ptr+two_isize_g2);
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000000000001111, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          const __m512 vt_i0_512 = _mm512_mul_ps(vM1_0000_1111_2222_3333,vf_i0_j0123);
          const __m512 vt_i1_512 = _mm512_mul_ps(vM1_0000_1111_2222_3333,vf_i1_j0123);
          const __m512 vt_i2_512 = _mm512_mul_ps(vM1_0000_1111_2222_3333,vf_i2_j0123);
          const __m512 vt_i3_512 = _mm512_mul_ps(vM1_0000_1111_2222_3333,vf_i3_j0123);

          __m512 vt0_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b00000000), vt_i0_512);
          __m512 vt1_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b01010101), vt_i1_512);
          __m512 vt2_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b10101010), vt_i2_512);
          __m512 vt3_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b11111111), vt_i3_512);



           //__m512 vt_512 = vt0_512;
           //vt0_512 = _mm512_add_ps(vt1_512, vt1_512);
           //vt2_512 = _mm512_add_ps(vt2_512, vt3_512);
           //vt_512 = _mm512_add_ps(vt0_512, vt2_512);
           //vt_512 = _mm512_mul_ps(vt_512, vM2_512);

            //__m512 vt_512 = vt0_512;
            //vt_512 = _mm512_add_ps(vt_512, vt1_512);
            //vt_512 = _mm512_add_ps(vt_512, vt2_512);
            //vt_512 = _mm512_add_ps(vt_512, vt3_512);
            //vt_512 = _mm512_mul_ps(vt_512, vM2_512);

           __m512 vt_512;
           vt0_512 = _mm512_add_ps(vt0_512, vt1_512);
           vt2_512 = _mm512_add_ps(vt2_512, vt3_512);
           vt_512 = _mm512_add_ps(vt0_512, vt2_512);
            vt_512 = _mm512_mul_ps(vt_512, vM2_512);

           //val[jj] = _mm512_reduce_add_ps (vt_512);
           query_values[i] = _mm512_reduce_add_ps (vt_512);
	  } //end jj loop
  //__m512 tmp = _mm512_loadu_ps(val);
  //_mm512_stream_ps (&query_values[ii*CHUNK],tmp);
	//} //end ii loop

	return;

}  // end of interp3_ghost_xyz_p

void ectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled) {

  const __m256  c1000 = _mm256_set_ps(-1.0,-0.0,-0.0,-0.0,-1.0,-0.0,-0.0,-0.0);
  const __m256  c2211 = _mm256_set_ps(-2.0,-2.0,-1.0,-1.0,-2.0,-2.0,-1.0,-1.0);
  const __m256  c3332 = _mm256_set_ps(-3.0,-3.0,-3.0,-2.0,-3.0,-3.0,-3.0,-2.0);
  const __m512  c1000_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-1.0,-0.0,-0.0,-0.0));
  const __m512  c2211_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-2.0,-2.0,-1.0,-1.0));
  const __m512  c3332_512 = _mm512_broadcast_f32x4(_mm_setr_ps(-3.0,-3.0,-3.0,-2.0));

  const __m512  c0010_512 = _mm512_set_ps(
      -0.0,-0.0,-0.0,-0.0,
      -0.0,-0.0,-0.0,-0.0,
      -1.0,-1.0,-1.0,-1.0,
      -0.0,-0.0,-0.0,-0.0
      );
  const __m512  c1122_512 = _mm512_set_ps(
      -1.0,-1.0,-1.0,-1.0,
      -1.0,-1.0,-1.0,-1.0,
      -2.0,-2.0,-2.0,-2.0,
      -2.0,-2.0,-2.0,-2.0
      );
  const __m512  c3233_512 = _mm512_set_ps(
      -3.0,-3.0,-3.0,-3.0,
      -2.0,-2.0,-2.0,-2.0,
      -3.0,-3.0,-3.0,-3.0,
      -3.0,-3.0,-3.0,-3.0
      );
  const __m512 vlagr_3333_2222_0000_1111_512  = _mm512_set_ps(
      -0.5,-0.5,-0.5,-0.5,
      +0.1666666667,0.1666666667,0.1666666667,0.1666666667,
      -0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,
       0.5,0.5,0.5,0.5
      );
  //print512(c1000_512,"512");
  //print512(c3332_512,"");
  //do{}while(1);

  const __m256 vlagr = _mm256_set_ps(-0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667);
  const __m512 vlagr_512 = _mm512_set_ps(
      -0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667,
      -0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667
      );
  const __m256  c33332222 = _mm256_set_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c22223333 = _mm256_setr_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c11110000 = _mm256_set_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  c00001111 = _mm256_setr_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  l0l1 = _mm256_set_ps (-0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,+0.5,+0.5,+0.5,+0.5);
  const __m256  l2l3 = _mm256_setr_ps(+0.1666666667,+0.1666666667,+0.1666666667,+0.1666666667,-0.5,-0.5,-0.5,-0.5);
  const int isize_g2 = isize_g[2];
  const int two_isize_g2 = 2*isize_g2;
  const int reg_plus = isize_g[1]*isize_g2 - two_isize_g2;
  const int NzNy = isize_g2 * isize_g[1];
  Real* Q_ptr = query_points;


  // std::cout << "KNL" << std::endl;
  //_mm_prefetch( (char*)Q_ptr,_MM_HINT_NTA);
#pragma omp parallel for
	for (int i = 0; i < N_pts; i++) {

#ifdef INTERP_USE_MORE_MEM_L1
		Real point[COORD_DIM];
		point[0] = Q_ptr[i*4+0];
		point[1] = Q_ptr[i*4+1];
		point[2] = Q_ptr[i*4+2];
    const int indxx = (int) Q_ptr[4*i + 3];
#else
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = Q_ptr[i*3+0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = Q_ptr[i*3+1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = Q_ptr[i*3+2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];
    // Q_ptr += 3;
		const int indxx = NzNy * grid_indx[0] + grid_indx[2] + isize_g2 * grid_indx[1] ;
#endif
    //_mm_prefetch( (char*)Q_ptr,_MM_HINT_T2);

    //__m256 vM0(vlagr), vM1(vlagr), vM2(vlagr);
    __m256 vM0(vlagr), vM2(vlagr);
    // __m256 vM0_tttt[4]; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest

    __m512 vM1_2222_3333_0000_1111(vlagr_3333_2222_0000_1111_512);
    __m512 vM2_512(vlagr_512);
    __m512 vM0_512(vlagr_512);

    {
    const __m256 vx0 = _mm256_set1_ps(point[0]);
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c1000));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c2211));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c3332));

    Real* dum, *dum2;
    const __m512 vx1_512 =  _mm512_set1_ps(point[1]);
    vM1_2222_3333_0000_1111  = _mm512_mul_ps(vM1_2222_3333_0000_1111 , _mm512_add_ps(vx1_512,c0010_512));
    vM1_2222_3333_0000_1111  = _mm512_mul_ps(vM1_2222_3333_0000_1111 , _mm512_add_ps(vx1_512,c1122_512));
    vM1_2222_3333_0000_1111  = _mm512_mul_ps(vM1_2222_3333_0000_1111 , _mm512_add_ps(vx1_512,c3233_512));

    const __m512 vx0_512 =  _mm512_set1_ps(point[0]);
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c1000_512));
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c2211_512));
    vM0_512  = _mm512_mul_ps(vM0_512 , _mm512_add_ps(vx0_512,c3332_512));

    const __m512 vx2_512 =  _mm512_set1_ps(point[2]);
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c1000_512));
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c2211_512));
    vM2_512  = _mm512_mul_ps(vM2_512 , _mm512_add_ps(vx2_512,c3332_512));

    const __m256 vx1 = _mm256_set1_ps(point[1]);
    __m256 tmp = _mm256_add_ps(vx1,c33332222); // x-3,...;x-2,...
    tmp = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c11110000));
    vM1_0000_1111 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c22223333));
    vM1_0000_1111 = _mm256_mul_ps(vM1_0000_1111, l0l1);
    vM1_2222_3333 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c00001111));
    vM1_2222_3333  = _mm256_mul_ps(vM1_2222_3333, l2l3);

    const __m256 vx2 = _mm256_set1_ps(point[2]);
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c1000));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c2211));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c3332));
    // todo remove permute completely by using different c's in the beginning
    vM2 = _mm256_permute_ps(vM2,0b00011011);
    ///print256(vM2,"256");
    ///print512(vM2_512,"256");
    ///print256(vM0,"256");
    ///print512(vM0_512,"256");
    ///do{}while(1);
    }

    int indx = 0;
    Real* reg_ptr = reg_grid_vals + indxx;//&reg_grid_vals[indxx];
    //_mm_prefetch( (char*)reg_ptr,_MM_HINT_T0);
		Real val = 0;

    // load all vfij
          __m256 vt;
          vt = _mm256_setzero_ps();

          __m512 vf_i0_j0123 = _mm512_setzero_ps();
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000000011110000, reg_ptr);
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000000000001111, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b1111000000000000, reg_ptr);
          vf_i0_j0123 = _mm512_mask_expandloadu_ps(vf_i0_j0123, 0b0000111100000000, reg_ptr+isize_g2);
           reg_ptr +=  reg_plus;

          __m512 vf_i1_j0123 = _mm512_setzero_ps();
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000000011110000, reg_ptr);
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000000000001111, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b1111000000000000, reg_ptr);
          vf_i1_j0123 = _mm512_mask_expandloadu_ps(vf_i1_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          __m512 vf_i2_j0123 = _mm512_setzero_ps();
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000000011110000, reg_ptr);
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000000000001111, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b1111000000000000, reg_ptr);
          vf_i2_j0123 = _mm512_mask_expandloadu_ps(vf_i2_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          __m512 vf_i3_j0123 = _mm512_setzero_ps();
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000000011110000, reg_ptr);
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000000000001111, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b1111000000000000, reg_ptr);
          vf_i3_j0123 = _mm512_mask_expandloadu_ps(vf_i3_j0123, 0b0000111100000000, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          reg_ptr = reg_grid_vals + indxx;//&reg_grid_vals[indxx];
          __m256 vf_i0_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          const __m256 vf_i0_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;
          const __m256 vf_i1_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          const __m256 vf_i1_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vf_i2_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          const __m256 vf_i2_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vf_i3_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          const __m256 vf_i3_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);


          const __m256 vt_i0_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i0_j01);
          const __m256 vt_i0_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i0_j23);
          const __m256 vt_i0 = _mm256_add_ps(vt_i0_j01, vt_i0_j23);
          const __m512 vt_i0_512 = _mm512_mul_ps(vM1_2222_3333_0000_1111,vf_i0_j0123);

          const __m256 vt_i1_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i1_j01);
          const __m256 vt_i1_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i1_j23);
          const __m256 vt_i1 = _mm256_add_ps(vt_i1_j01, vt_i1_j23);
          const __m512 vt_i1_512 = _mm512_mul_ps(vM1_2222_3333_0000_1111,vf_i1_j0123);

          const __m256 vt_i2_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i2_j01);
          const __m256 vt_i2_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i2_j23);
          const __m256 vt_i2 = _mm256_add_ps(vt_i2_j01, vt_i2_j23);
          const __m512 vt_i2_512 = _mm512_mul_ps(vM1_2222_3333_0000_1111,vf_i2_j0123);

          const __m256 vt_i3_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i3_j01);
          const __m256 vt_i3_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i3_j23);
          const __m256 vt_i3 = _mm256_add_ps(vt_i3_j01, vt_i3_j23);
          const __m512 vt_i3_512 = _mm512_mul_ps(vM1_2222_3333_0000_1111,vf_i3_j0123);

          const __m256 vt0 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b11111111), vt_i0);
          const __m256 vt1 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b10101010), vt_i1);
          const __m256 vt2 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b01010101), vt_i2);
          const __m256 vt3 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b00000000), vt_i3);

          const __m512 vt0_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b00000000), vt_i0_512);
          const __m512 vt1_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b01010101), vt_i1_512);
          const __m512 vt2_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b10101010), vt_i2_512);
          const __m512 vt3_512 = _mm512_mul_ps(_mm512_permute_ps(vM0_512,0b11111111), vt_i3_512);


           __m512 vt_512 = vt0_512;
           vt_512 = _mm512_add_ps(vt_512, vt1_512);
           vt_512 = _mm512_add_ps(vt_512, vt2_512);
           vt_512 = _mm512_add_ps(vt_512, vt3_512);
           vt_512 = _mm512_mul_ps(vt_512, vM2_512);


           vt = _mm256_add_ps(vt0, vt1);
           vt = _mm256_add_ps(vt, vt2);
           vt = _mm256_add_ps(vt, vt3);

           vt = _mm256_mul_ps(vM2, vt);
            val = sum8(vt);
		        query_values[i] = val;

           //print256(vt,"");
           //print512(vt_512,"");
           //do{}while(1);
           val = _mm512_reduce_add_ps (vt_512);
		       query_values[i] = val;
	}

	return;

}  // end of interp3_ghost_xyz_p

#elif defined(__AVX2__) || defined(HASWELL)
void vectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled) {

#ifdef INTERP_DEBUG
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  PCOUT << "In Haswel kernel\n";
#endif
  const __m256  c1000 = _mm256_set_ps(-1.0,-0.0,-0.0,-0.0,-1.0,-0.0,-0.0,-0.0);
  const __m256  c2211 = _mm256_set_ps(-2.0,-2.0,-1.0,-1.0,-2.0,-2.0,-1.0,-1.0);
  const __m256  c3332 = _mm256_set_ps(-3.0,-3.0,-3.0,-2.0,-3.0,-3.0,-3.0,-2.0);

  const __m256 vlagr = _mm256_set_ps(-0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667);
  const __m256  c33332222 = _mm256_set_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c22223333 = _mm256_setr_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c11110000 = _mm256_set_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  c00001111 = _mm256_setr_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  l0l1 = _mm256_set_ps (-0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,+0.5,+0.5,+0.5,+0.5);
  const __m256  l2l3 = _mm256_setr_ps(+0.1666666667,+0.1666666667,+0.1666666667,+0.1666666667,-0.5,-0.5,-0.5,-0.5);
  const int isize_g2 = isize_g[2];
  const int two_isize_g2 = 2*isize_g2;
  const int three_isize_g2 = 3*isize_g2;
  const int reg_plus = isize_g[1]*isize_g2;
  const int NzNy = isize_g2 * isize_g[1];
  Real* Q_ptr = query_points;
  //std::cout << "AVX2" << std::endl;
  //_mm_prefetch( (char*)Q_ptr,_MM_HINT_NTA);

 #pragma omp parallel for
	 for (int i = 0; i < N_pts; i++) {
//  int CHUNK=8;
//#pragma omp parallel for
//  for (int ii = 0; ii < (int)std::ceil(N_pts/(float)CHUNK); ii++) {
//		Real val[CHUNK];
//	for (int jj = 0; jj < CHUNK; jj++) {
//    int i = ii*CHUNK + jj;
#ifdef INTERP_USE_MORE_MEM_L1
		Real point[COORD_DIM];
		point[0] = Q_ptr[i*4+0];
		point[1] = Q_ptr[i*4+1];
		point[2] = Q_ptr[i*4+2];
    const int indxx = (int) Q_ptr[4*i + 3];
#else
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = Q_ptr[i*3+0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = Q_ptr[i*3+1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = Q_ptr[i*3+2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];
    // Q_ptr += 3;
		const int indxx = NzNy * grid_indx[0] + grid_indx[2] + isize_g2 * grid_indx[1] ;
#endif

    ////_mm_prefetch( (char*)Q_ptr,_MM_HINT_T2);
//    int indx = 0;
    Real* reg_ptr = reg_grid_vals + indxx;//&reg_grid_vals[indxx];
    //_mm_prefetch( (char*)reg_ptr,_MM_HINT_T0);



    //__m256 vM0(vlagr), vM1(vlagr), vM2(vlagr);
    __m256 vM0(vlagr), vM2(vlagr);
    // __m256 vM0_tttt[4]; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest


    {
    const __m256 vx0 = _mm256_set1_ps(point[0]);
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c1000));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c2211));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c3332));

    const __m256 vx1 = _mm256_set1_ps(point[1]);
    __m256 tmp = _mm256_add_ps(vx1,c33332222); // x-3,...;x-2,...
    tmp = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c11110000));
    vM1_0000_1111 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c22223333));
    vM1_0000_1111 = _mm256_mul_ps(vM1_0000_1111, l0l1);
    vM1_2222_3333 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c00001111));
    vM1_2222_3333  = _mm256_mul_ps(vM1_2222_3333, l2l3);
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c1000));
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c2211));
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c3332));

    const __m256 vx2 = _mm256_set1_ps(point[2]);
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c1000));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c2211));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c3332));
    // todo remove permute completely by using different c's in the beginning
    vM2 = _mm256_permute_ps(vM2,0b00011011);
    //vM2 = _mm256_shuffle_ps(vM2,vM2,_MM_SHUFFLE(0, 1, 2, 3));

    //const Real* M1 = (Real*)&vM1;
    //vM1_0000_1111 = _mm256_set_ps(M1[7],M1[7],M1[7],M1[7],M1[6],M1[6],M1[6],M1[6]);
    //vM1_2222_3333 = _mm256_set_ps(M1[5],M1[5],M1[5],M1[5],M1[4],M1[4],M1[4],M1[4]);
    // vM1_0000_1111 = _mm256_set_ps(M1[3],M1[3],M1[3],M1[3],M1[2],M1[2],M1[2],M1[2]);
    // vM1_2222_3333 = _mm256_set_ps(M1[1],M1[1],M1[1],M1[1],M1[0],M1[0],M1[0],M1[0]);
    //vM0_tttt[0] = _mm256_permute_ps(vM0,0b11111111); // last element
    //vM0_tttt[1] = _mm256_permute_ps(vM0,0b10101010);
    //vM0_tttt[2] = _mm256_permute_ps(vM0,0b01010101);
    //vM0_tttt[3] = _mm256_permute_ps(vM0,0b00000000);
    }


    // load all vfij
          __m256 vt;
          vt = _mm256_setzero_ps();
          const __m256 vf_i0_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          const __m256 vf_i0_j23 = _mm256_loadu2_m128(reg_ptr+two_isize_g2, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vf_i1_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          const __m256 vf_i1_j23 = _mm256_loadu2_m128(reg_ptr+two_isize_g2, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vf_i2_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          const __m256 vf_i2_j23 = _mm256_loadu2_m128(reg_ptr+two_isize_g2, reg_ptr+three_isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vf_i3_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          const __m256 vf_i3_j23 = _mm256_loadu2_m128(reg_ptr+two_isize_g2, reg_ptr+three_isize_g2);

          //__m256 vt0, vt1, vt2, vt3;
          //vt0 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b11111111),_mm256_fmadd_ps(vM1_0000_1111, vf_i0_j01, _mm256_mul_ps(vM1_2222_3333, vf_i0_j23)));
          //vt1 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b10101010),_mm256_fmadd_ps(vM1_0000_1111, vf_i1_j01, _mm256_mul_ps(vM1_2222_3333, vf_i1_j23)));
          //vt2 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b01010101),_mm256_fmadd_ps(vM1_0000_1111, vf_i2_j01, _mm256_mul_ps(vM1_2222_3333, vf_i2_j23)));
          //vt3 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b00000000),_mm256_fmadd_ps(vM1_0000_1111, vf_i3_j01, _mm256_mul_ps(vM1_2222_3333, vf_i3_j23)));

          vt = _mm256_fmadd_ps(_mm256_permute_ps(vM0,0b11111111),_mm256_fmadd_ps(vM1_0000_1111, vf_i0_j01, _mm256_mul_ps(vM1_2222_3333, vf_i0_j23)) , vt);
          vt = _mm256_fmadd_ps(_mm256_permute_ps(vM0,0b10101010),_mm256_fmadd_ps(vM1_0000_1111, vf_i1_j01, _mm256_mul_ps(vM1_2222_3333, vf_i1_j23)) , vt);
          vt = _mm256_fmadd_ps(_mm256_permute_ps(vM0,0b01010101),_mm256_fmadd_ps(vM1_0000_1111, vf_i2_j01, _mm256_mul_ps(vM1_2222_3333, vf_i2_j23)) , vt);
          vt = _mm256_fmadd_ps(_mm256_permute_ps(vM0,0b00000000),_mm256_fmadd_ps(vM1_0000_1111, vf_i3_j01, _mm256_mul_ps(vM1_2222_3333, vf_i3_j23)) , vt);

          vt = _mm256_mul_ps(vM2, vt);
          //val[jj] = sum8(vt);
          query_values[i] = sum8(vt);
	        }

//  __m256 tmp = _mm256_loadu_ps(val);
//  _mm256_stream_ps (&query_values[ii*CHUNK],tmp);
//	}

	return;

}  // end of interp3_ghost_xyz_p

#else
void vectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled) {

#ifdef INTERP_DEBUG
	int nprocs, procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  PCOUT << "In ivybridge  kernel\n";
#endif
  const __m256  c1000 = _mm256_set_ps(-1.0,-0.0,-0.0,-0.0,-1.0,-0.0,-0.0,-0.0);
  const __m256  c2211 = _mm256_set_ps(-2.0,-2.0,-1.0,-1.0,-2.0,-2.0,-1.0,-1.0);
  const __m256  c3332 = _mm256_set_ps(-3.0,-3.0,-3.0,-2.0,-3.0,-3.0,-3.0,-2.0);

  const __m256 vlagr = _mm256_set_ps(-0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667);
  const __m256  c33332222 = _mm256_set_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c22223333 = _mm256_setr_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
  const __m256  c11110000 = _mm256_set_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  c00001111 = _mm256_setr_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
  const __m256  l0l1 = _mm256_set_ps (-0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,+0.5,+0.5,+0.5,+0.5);
  const __m256  l2l3 = _mm256_setr_ps(+0.1666666667,+0.1666666667,+0.1666666667,+0.1666666667,-0.5,-0.5,-0.5,-0.5);
  const int isize_g2 = isize_g[2];
  const int two_isize_g2 = 2*isize_g2;
  const int reg_plus = isize_g[1]*isize_g2 - two_isize_g2;
  const int NzNy = isize_g2 * isize_g[1];
  Real* Q_ptr = query_points;
  //_mm_prefetch( (char*)Q_ptr,_MM_HINT_NTA);
	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = Q_ptr[0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = Q_ptr[1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = Q_ptr[2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];
    Q_ptr += 3;

		const int indxx = NzNy * grid_indx[0] + grid_indx[2] + isize_g2 * grid_indx[1] ;
    //_mm_prefetch( (char*)Q_ptr,_MM_HINT_T2);

//    int indx = 0;
    Real* reg_ptr = reg_grid_vals + indxx;//&reg_grid_vals[indxx];
    //_mm_prefetch( (char*)reg_ptr,_MM_HINT_T0);
		Real val = 0;



//  _mm_prefetch( (char*)reg_ptr,_MM_HINT_T1 );
//  _mm_prefetch( (char*)reg_ptr+isize_g2,_MM_HINT_T1 );
//  reg_ptr += two_isize_g2;
//  _mm_prefetch( (char*)reg_ptr,_MM_HINT_T0 );
//  _mm_prefetch( (char*)reg_ptr+isize_g2,_MM_HINT_T0 );
//  reg_ptr += reg_plus;
//           const __m256 vf_i0_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i0_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i1_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i1_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i2_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i2_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i3_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i3_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           // reg_ptr +=  reg_plus;




    //__m256 vM0(vlagr), vM1(vlagr), vM2(vlagr);
    __m256 vM0(vlagr), vM2(vlagr);
    // __m256 vM0_tttt[4]; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest
    __m256 vM0_tttt[4];


    {
    const __m256 vx0 = _mm256_set1_ps(point[0]);
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c1000));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c2211));
    vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c3332));

    const __m256 vx1 = _mm256_set1_ps(point[1]);
    __m256 tmp = _mm256_add_ps(vx1,c33332222); // x-3,...;x-2,...
    tmp = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c11110000));
    vM1_0000_1111 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c22223333));
    vM1_0000_1111 = _mm256_mul_ps(vM1_0000_1111, l0l1);
    vM1_2222_3333 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c00001111));
    vM1_2222_3333  = _mm256_mul_ps(vM1_2222_3333, l2l3);
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c1000));
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c2211));
    //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c3332));

    const __m256 vx2 = _mm256_set1_ps(point[2]);
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c1000));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c2211));
    vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c3332));
    // todo remove permute completely by using different c's in the beginning
    vM2 = _mm256_permute_ps(vM2,0b00011011);
    //vM2 = _mm256_shuffle_ps(vM2,vM2,_MM_SHUFFLE(0, 1, 2, 3));

    //const Real* M1 = (Real*)&vM1;
    //vM1_0000_1111 = _mm256_set_ps(M1[7],M1[7],M1[7],M1[7],M1[6],M1[6],M1[6],M1[6]);
    //vM1_2222_3333 = _mm256_set_ps(M1[5],M1[5],M1[5],M1[5],M1[4],M1[4],M1[4],M1[4]);
    // vM1_0000_1111 = _mm256_set_ps(M1[3],M1[3],M1[3],M1[3],M1[2],M1[2],M1[2],M1[2]);
    // vM1_2222_3333 = _mm256_set_ps(M1[1],M1[1],M1[1],M1[1],M1[0],M1[0],M1[0],M1[0]);
    vM0_tttt[0] = _mm256_permute_ps(vM0,0b11111111); // last element
    vM0_tttt[1] = _mm256_permute_ps(vM0,0b10101010);
    vM0_tttt[2] = _mm256_permute_ps(vM0,0b01010101);
    vM0_tttt[3] = _mm256_permute_ps(vM0,0b00000000);
    }


    // load all vfij




          //
          const __m256 vf_i0_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;

          const __m256 vt_i0_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i0_j01); //8

          const __m256 vf_i0_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;
          const __m256 vt_i0_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i0_j23);//8

          const __m256 vt_i0 = _mm256_add_ps(vt_i0_j01, vt_i0_j23);//8

          //
          const __m256 vf_i1_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;

          const __m256 vt_i1_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i1_j01);//8

          const __m256 vf_i1_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;
          const __m256 vt_i1_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i1_j23);//8

          const __m256 vt_i1 = _mm256_add_ps(vt_i1_j01, vt_i1_j23);//8

          //
          const __m256 vf_i2_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;

          const __m256 vt_i2_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i2_j01);//8

          const __m256 vf_i2_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr +=  reg_plus;

          const __m256 vt_i2_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i2_j23);//8
          const __m256 vt_i2 = _mm256_add_ps(vt_i2_j01, vt_i2_j23);//8

          //
          const __m256 vf_i3_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          // reg_ptr +=  reg_plus;

          const __m256 vt_i3_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i3_j01);//8
          const __m256 vf_i3_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          const __m256 vt_i3_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i3_j23);//8
          const __m256 vt_i3 = _mm256_add_ps(vt_i3_j01, vt_i3_j23);//8

          const __m256 vt0 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b11111111), vt_i0);//8
          const __m256 vt1 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b10101010), vt_i1);//8
          const __m256 vt2 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b01010101), vt_i2);//8
          const __m256 vt3 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b00000000), vt_i3);//8
          //const __m256 vt0 = _mm256_mul_ps(vM0_tttt[0], vt_i0);
          //const __m256 vt1 = _mm256_mul_ps(vM0_tttt[1], vt_i1);
          //const __m256 vt2 = _mm256_mul_ps(vM0_tttt[2], vt_i2);
          //const __m256 vt3 = _mm256_mul_ps(vM0_tttt[3], vt_i3);

          __m256 vt = _mm256_add_ps(vt0, vt1);//8
          vt = _mm256_add_ps(vt, vt2);//8
          vt = _mm256_add_ps(vt, vt3);//8

          vt = _mm256_mul_ps(vM2, vt);//8
          val = sum8(vt);//7
		      query_values[i] = val;
	}

	return;

}  // end of interp3_ghost_xyz_p
#endif

#endif



// void vectorized_interp3_ghost_xyz_p(Real* __restrict reg_grid_vals, int data_dof, const int* __restrict N_reg,
// 		const int* __restrict N_reg_g, const int * __restrict isize_g, const int* __restrict istart, const int N_pts,
// 		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
// 		bool query_values_already_scaled) {
//
//   const __m256  c1000 = _mm256_set_ps(-1.0,-0.0,-0.0,-0.0,-1.0,-0.0,-0.0,-0.0);
//   const __m256  c2211 = _mm256_set_ps(-2.0,-2.0,-1.0,-1.0,-2.0,-2.0,-1.0,-1.0);
//   const __m256  c3332 = _mm256_set_ps(-3.0,-3.0,-3.0,-2.0,-3.0,-3.0,-3.0,-2.0);
//
//   const __m256 vlagr = _mm256_set_ps(-0.1666666667,0.5,-0.5, 0.1666666667,-0.1666666667,0.5,-0.5, 0.1666666667);
//   const __m256  c33332222 = _mm256_set_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
//   const __m256  c22223333 = _mm256_setr_ps(-3.0,-3.0,-3.0,-3.0,-2.0,-2.0,-2.0,-2.0);
//   const __m256  c11110000 = _mm256_set_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
//   const __m256  c00001111 = _mm256_setr_ps(-1.0,-1.0,-1.0,-1.0,0,0,0,0);
//   const __m256  l0l1 = _mm256_set_ps (-0.1666666667,-0.1666666667,-0.1666666667,-0.1666666667,+0.5,+0.5,+0.5,+0.5);
//   const __m256  l2l3 = _mm256_setr_ps(+0.1666666667,+0.1666666667,+0.1666666667,+0.1666666667,-0.5,-0.5,-0.5,-0.5);
// 	for (int i = 0; i < N_pts; i++) {
// 		Real point[COORD_DIM];
// 		int grid_indx[COORD_DIM];
//
// 		point[0] = query_points[COORD_DIM * i + 0] * N_reg_g[0];
// 		grid_indx[0] = ((int)(point[0])) - 1;
// 		point[0] -= grid_indx[0];
//
// 		point[1] = query_points[COORD_DIM * i + 1] * N_reg_g[1];
// 		grid_indx[1] = ((int)(point[1])) - 1;
// 		point[1] -= grid_indx[1];
//
// 		point[2] = query_points[COORD_DIM * i + 2] * N_reg_g[2];
// 		grid_indx[2] = ((int)(point[2])) - 1;
// 		point[2] -= grid_indx[2];
//
// 		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
//     Real* reg_ptr = &reg_grid_vals[indxx];
// 		Real val = 0;
//     int indx = 0;
//     const int isize_g2 = isize_g[2];
//     const int two_isize_g2 = 2*isize_g[2];
//     const int reg_plus = isize_g[1]*isize_g2 - two_isize_g2;
//
//
//     __m256 vM0(vlagr), vM1(vlagr), vM2(vlagr);
//     // __m256 vM0_tttt[4]; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
//     __m256 vM1_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
//     __m256 vM1_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest
//     __m256 vM0_tttt[4];
//
//
//     {
//     const __m256 vx0 = _mm256_set1_ps(point[0]);
//     const __m256 vx1 = _mm256_set1_ps(point[1]);
//     const __m256 vx2 = _mm256_set1_ps(point[2]);
//     vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c1000));
//     vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c2211));
//     vM0 = _mm256_mul_ps(vM0, _mm256_add_ps(vx0,c3332));
//
//     __m256 tmp = _mm256_add_ps(vx1,c33332222); // x-3,...;x-2,...
//     tmp = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c11110000));
//     vM1_0000_1111 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c22223333));
//     vM1_0000_1111 = _mm256_mul_ps(vM1_0000_1111, l0l1);
//     vM1_2222_3333 = _mm256_mul_ps(tmp, _mm256_add_ps(vx1,c00001111));
//     vM1_2222_3333  = _mm256_mul_ps(vM1_2222_3333, l2l3);
//     //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c1000));
//     //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c2211));
//     //vM1 = _mm256_mul_ps(vM1, _mm256_add_ps(vx1,c3332));
//
//     vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c1000));
//     vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c2211));
//     vM2 = _mm256_mul_ps(vM2, _mm256_add_ps(vx2,c3332));
//     // todo remove permute completely by using different c's in the beginning
//     vM2 = _mm256_permute_ps(vM2,0b00011011);
//     //vM2 = _mm256_shuffle_ps(vM2,vM2,_MM_SHUFFLE(0, 1, 2, 3));
//
//     //const Real* M1 = (Real*)&vM1;
//     //vM1_0000_1111 = _mm256_set_ps(M1[7],M1[7],M1[7],M1[7],M1[6],M1[6],M1[6],M1[6]);
//     //vM1_2222_3333 = _mm256_set_ps(M1[5],M1[5],M1[5],M1[5],M1[4],M1[4],M1[4],M1[4]);
//     // vM1_0000_1111 = _mm256_set_ps(M1[3],M1[3],M1[3],M1[3],M1[2],M1[2],M1[2],M1[2]);
//     // vM1_2222_3333 = _mm256_set_ps(M1[1],M1[1],M1[1],M1[1],M1[0],M1[0],M1[0],M1[0]);
//     vM0_tttt[0] = _mm256_permute_ps(vM0,0b11111111); // last element
//     vM0_tttt[1] = _mm256_permute_ps(vM0,0b10101010);
//     vM0_tttt[2] = _mm256_permute_ps(vM0,0b01010101);
//     vM0_tttt[3] = _mm256_permute_ps(vM0,0b00000000);
//     }
//
//
//     // load all vfij
//           const __m256 vf_i0_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i0_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i1_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i1_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i2_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i2_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr +=  reg_plus;
//
//           const __m256 vf_i3_j01 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           reg_ptr += two_isize_g2;
//           const __m256 vf_i3_j23 = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
//           // reg_ptr +=  reg_plus;
//
//           const __m256 vt_i0_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i0_j01);
//           const __m256 vt_i0_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i0_j23);
//           const __m256 vt_i0 = _mm256_add_ps(vt_i0_j01, vt_i0_j23);
//
//           const __m256 vt_i1_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i1_j01);
//           const __m256 vt_i1_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i1_j23);
//           const __m256 vt_i1 = _mm256_add_ps(vt_i1_j01, vt_i1_j23);
//
//           const __m256 vt_i2_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i2_j01);
//           const __m256 vt_i2_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i2_j23);
//           const __m256 vt_i2 = _mm256_add_ps(vt_i2_j01, vt_i2_j23);
//
//           const __m256 vt_i3_j01 = _mm256_mul_ps(vM1_0000_1111, vf_i3_j01);
//           const __m256 vt_i3_j23 = _mm256_mul_ps(vM1_2222_3333, vf_i3_j23);
//           const __m256 vt_i3 = _mm256_add_ps(vt_i3_j01, vt_i3_j23);
//
//           const __m256 vt0 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b11111111), vt_i0);
//           const __m256 vt1 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b10101010), vt_i1);
//           const __m256 vt2 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b01010101), vt_i2);
//           const __m256 vt3 = _mm256_mul_ps(_mm256_permute_ps(vM0,0b00000000), vt_i3);
//           //const __m256 vt0 = _mm256_mul_ps(vM0_tttt[0], vt_i0);
//           //const __m256 vt1 = _mm256_mul_ps(vM0_tttt[1], vt_i1);
//           //const __m256 vt2 = _mm256_mul_ps(vM0_tttt[2], vt_i2);
//           //const __m256 vt3 = _mm256_mul_ps(vM0_tttt[3], vt_i3);
//
//           __m256 vt = _mm256_add_ps(vt0, vt1);
//           vt = _mm256_add_ps(vt, vt2);
//           vt = _mm256_add_ps(vt, vt3);
//
//           vt = _mm256_mul_ps(vM2, vt);
//           val = sum8(vt);
// 		      query_values[i] = val;
// 	}
//
// 	return;
//
// }  // end of interp3_ghost_xyz_p


void rescale_xyz(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		const int N_pts, Real* Q_) {

	if (g_size == 0)
		return;
	Real hp[3];
	Real h[3];
	hp[0] = 1. / N_reg_g[0]; // New mesh size
	hp[1] = 1. / N_reg_g[1]; // New mesh size
	hp[2] = 1. / N_reg_g[2]; // New mesh size

	h[0] = 1. / (N_reg[0]); // old mesh size
	h[1] = 1. / (N_reg[1]); // old mesh size
	h[2] = 1. / (N_reg[2]); // old mesh size

	const Real factor0 = (1. - (2. * g_size + 1.) * hp[0]) / (1. - h[0]);
	const Real factor1 = (1. - (2. * g_size + 1.) * hp[1]) / (1. - h[1]);
	const Real factor2 = (1. - (2. * g_size + 1.) * hp[2]) / (1. - h[2]);
  const Real iX0 = istart[0]*h[0];
  const Real iX1 = istart[1]*h[1];
  const Real iX2 = istart[2]*h[2];

	for (int i = 0; i < N_pts; i++) {
		Q_[0 + COORD_DIM * i] = (Q_[0 + COORD_DIM * i]
				- iX0) * factor0 + g_size * hp[0];
		Q_[1 + COORD_DIM * i] = (Q_[1 + COORD_DIM * i]
				- iX1) * factor1 + g_size * hp[1];
		Q_[2 + COORD_DIM * i] = (Q_[2 + COORD_DIM * i]
				- iX2) * factor2 + g_size * hp[2];
	}
	return;
} // end of rescale_xyz


#ifdef FAST_INTERPV

//#include "v1.cpp" // corresponding optimized version
void vec_torized_interp3_ghost_xyz_p( Real* __restrict reg_grid_vals, int data_dof, const int* N_reg,
		const int* N_reg_g, const int * isize_g, const int* istart, const int N_pts,
		const int g_size, Real* __restrict query_points, Real* __restrict query_values,
		bool query_values_already_scaled) {

  const __m128  c1000 = _mm_set_ps(-1.0,-0.0,-0.0,-0.0);
  const __m128  c2211 = _mm_set_ps(-2.0,-2.0,-1.0,-1.0);
  const __m128  c3332 = _mm_set_ps(-3.0,-3.0,-3.0,-2.0);
  const __m128 vlagr = _mm_set_ps(-0.1666666667,0.5,-0.5, 0.1666666667);

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = query_points[COORD_DIM * i + 0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = query_points[COORD_DIM * i + 1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = query_points[COORD_DIM * i + 2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];

		// for (int j = 0; j < COORD_DIM; j++) {
		//	Real x = point[j];
		//	for (int k = 0; k < 4; k++) {
		//		M[j][k] = lagr_denom[k];
		//		for (int l = 0; l < 4; l++) {
		//			if (k != l)
		//				M[j][k] *= (x - l);
		//		}
		//	}
		// }
    //M[0][0] = lagr_denom[0];
    //M[0][1] = lagr_denom[1];
    //M[0][2] = lagr_denom[2];
    //M[0][3] = lagr_denom[3];

    //M[0][0] *= (x-1);
    //M[0][1] *= (x-0);
    //M[0][2] *= (x-0);
    //M[0][3] *= (x-0);

    //M[0][0] *= (x-2);
    //M[0][1] *= (x-2);
    //M[0][2] *= (x-1);
    //M[0][3] *= (x-1);

    //M[0][0] *= (x-3);
    //M[0][1] *= (x-3);
    //M[0][2] *= (x-3);
    //M[0][3] *= (x-2);

    __m128 vx;

    __m128 vM0(vlagr), vM1(vlagr), vM2(vlagr);
    __m256 vM0_tttt[4]; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM1_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest
          __m256 vVal_, vPtr_;
          __m256 vVal2_;//= _mm256_set1_ps(point[0]);


    vx = _mm_set1_ps(point[0]);
    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c1000));
    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c2211));
    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c3332));

    vx = _mm_set1_ps(point[1]);
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c1000));
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c2211));
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c3332));

    vx = _mm_set1_ps(point[2]);
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c1000));
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c2211));
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c3332));
    vM2 = _mm_shuffle_ps(vM2,vM2,_MM_SHUFFLE(0, 1, 2, 3));


    vM1_0000_1111 = _mm256_set_m128(
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(3, 3, 3, 3)), // M[1][0]
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(2, 2, 2, 2)));// M[1][1]
    vM1_2222_3333 = _mm256_set_m128(
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(1, 1, 1, 1)), // M[1][2]
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(0, 0, 0, 0)));// M[1][3]


    Real* M0 = (Real*)&vM0;
    vM0_tttt[3] = _mm256_set1_ps(M0[0]);
    vM0_tttt[2] = _mm256_set1_ps(M0[1]);
    vM0_tttt[1] = _mm256_set1_ps(M0[2]);
    vM0_tttt[0] = _mm256_set1_ps(M0[3]);

    //vM0_tttt[0] = _mm256_set_m128(
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(3, 3, 3, 3)), // M[0][0]
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(3, 3, 3, 3)));// M[0][0]
    //vM0_tttt[1] = _mm256_set_m128(
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(2, 2, 2, 2)), // M[0][0]
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(2, 2, 2, 2)));// M[0][0]
    //vM0_tttt[2] = _mm256_set_m128(
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(1, 1, 1, 1)), // M[0][0]
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(1, 1, 1, 1)));// M[0][0]
    //vM0_tttt[3] = _mm256_set_m128(
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(0, 0, 0, 0)), // M[0][0]
    //          _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(0, 0, 0, 0)));// M[0][0]
    //Real* dum1 = (Real*)&vM0;
    //Real* dum2 = (Real*)&vM1;
    //Real* dum3 = (Real*)&vM2;
    //Real* dum4 = (Real*)&vM1_0000_1111;
    //Real* dum5 = (Real*)&vM1_2222_3333;
    //Real* dum6 = (Real*)&vVal_;
    //Real* dum7 = (Real*)&vM0_tttt[0];
    //Real* dum8 = (Real*)&vM0_tttt[1];
    //Real* dum9 = (Real*)&vM0_tttt[2];
    //Real* dum10 = (Real*)&vM0_tttt[3];
    //query_values[i] = point[0]*point[1]*point[2];
    //continue;
		//query_values[i] = dum1[0]*dum2[0]*dum3[2]*dum4[5]*dum5[5]*dum5[0]*dum5[2]
    //  *dum7[0]*dum8[5]*dum9[0]*dum10[2];//*dum5[3];
    //continue;


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
    Real* reg_ptr = &reg_grid_vals[indxx];
		Real val = 0;
    //int indx = 0;
    const int isize_g2 = isize_g[2];
    const int two_isize_g2 = 2*isize_g[2];
    const int reg_plus = isize_g[1]*isize_g2 - two_isize_g2;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          vVal_ = _mm256_setzero_ps();


          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM1_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[0]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM1_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[0]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);


          reg_ptr +=  reg_plus;

      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM1_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[1]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);


          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM1_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[1]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);
          reg_ptr +=  reg_plus;

      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM1_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[2]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM1_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[2]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          reg_ptr +=  reg_plus;
      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM1_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[3]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM1_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_tttt[3]);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vm_inv = M[2][0], [1] [2] [3] in reverse order
          __m256 vM2_256 = _mm256_set_m128(vM2, vM2);
          vVal_ = _mm256_mul_ps(vVal_, vM2_256);
          val = sum8(vVal_);
		query_values[i] = val;
	}

	return;

}  // end of interp3_ghost_xyz_p


void _vectorized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
    //__m256 vM_ = _mm256_set_ps(M[2][0], M[2][1], M[2][2], M[2][3], 0, 0, 0, 0);
    register __m128 vM_ = _mm_set_ps(M[2][0], M[2][1], M[2][2], M[2][3]);
    //__m128 vM_ = _mm_loadu_ps(&M[2][0]);
    register __m128 vVal =  _mm_setzero_ps();
    register __m128 vVal_;
    // std::cout << "indxx = " << indxx << std::endl;
//		Real val = 0;
    int indx = 0;
		for (int j0 = 0; j0 < 4; j0++) {
			for (int j1 = 0; j1 < 4; j1++) {
          const register __m128 M0M1 = _mm_set1_ps(M[0][j0]*M[1][j1]);

          __m128 vPtr_ = _mm_loadu_ps(&reg_grid_vals[indx + indxx]);
          vVal_ = _mm_mul_ps(vM_, vPtr_);
          vVal_ = _mm_mul_ps(vVal_, M0M1);

          //__m256 vVal_;
          //__m256 vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], 0, 0, 0, 0);
          //vVal_ = _mm256_mul_ps(vM_, vPtr_);
          //Real val_ = sum8(vVal_);
          //val += (val_[4] + val_[5] + val_[6] + val_[7]) * M0M1;
          vVal = _mm_add_ps(vVal, vVal_);
          indx += isize_g[2];
			}
      indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];
		}
    vVal = _mm_hadd_ps(vVal, vVal);
    vVal = _mm_hadd_ps(vVal, vVal);
    Real* val_ = (Real*)&vVal;
		query_values[i] = val_[0];
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void __vectorized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
//		Real val = 0;
    int indx = 0;


		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1, M0M1_2;
          __m256 vVal_, vM0M1_, vM2_, vM_, vPtr_;
          Real* ptr, *ptr2;

					// val_ = M[2][0] * ptr[0];
					// val_ += M[2][1] * ptr[1];
					// val_ += M[2][2] * ptr[2];
					// val_ += M[2][3] * ptr[3];
          // val += val_ * M0M1;
          // indx += isize_g[2];
          // M0M1 = M[0][0]*M[1][1];
          // ptr = &reg_grid_vals[indx + indxx];
					// val_ = M[2][0] * ptr[0];
					// val_ += M[2][1] * ptr[1];
					// val_ += M[2][2] * ptr[2];
					// val_ += M[2][3] * ptr[3];
          // val += val_ * M0M1;
          //indx += isize_g[2];

          M0M1 = M[0][0]*M[1][0];
          M0M1_2 = M[0][0]*M[1][1];
          vM2_ = _mm256_set_ps(M[2][0], M[2][1], M[2][2], M[2][3], M[2][0], M[2][1], M[2][2], M[2][3]);
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);
          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          //----//

          M0M1 = M[0][0]*M[1][2];
          M0M1_2 = M[0][0]*M[1][3];
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          indx += isize_g[2];
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);

          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];


          // M0M1 = M[0][0]*M[1][2];
          // ptr = &reg_grid_vals[indx + indxx];
					// val_ = M[2][0] * ptr[0];
					// val_ += M[2][1] * ptr[1];
					// val_ += M[2][2] * ptr[2];
					// val_ += M[2][3] * ptr[3];
          // val += val_ * M0M1;
          // indx += isize_g[2];

          // M0M1 = M[0][0]*M[1][3];
          // ptr = &reg_grid_vals[indx + indxx];
					// val_ = M[2][0] * ptr[0];
					// val_ += M[2][1] * ptr[1];
					// val_ += M[2][2] * ptr[2];
					// val_ += M[2][3] * ptr[3];
          // val += val_ * M0M1;
          // indx += isize_g[2];
          // indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][1]*M[1][0];
          M0M1_2 = M[0][1]*M[1][1];
          vM2_ = _mm256_set_ps(M[2][0], M[2][1], M[2][2], M[2][3], M[2][0], M[2][1], M[2][2], M[2][3]);
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);
          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          //----//

          M0M1 = M[0][1]*M[1][2];
          M0M1_2 = M[0][1]*M[1][3];
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          indx += isize_g[2];
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);

          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

          //M0M1 = M[0][1]*M[1][0];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][1]*M[1][1];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][1]*M[1][2];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];

          //M0M1 = M[0][1]*M[1][3];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];
          //indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // -//----------------------------------- //
          M0M1 = M[0][2]*M[1][0];
          M0M1_2 = M[0][2]*M[1][1];
          vM2_ = _mm256_set_ps(M[2][0], M[2][1], M[2][2], M[2][3], M[2][0], M[2][1], M[2][2], M[2][3]);
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);
          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          //----//

          M0M1 = M[0][2]*M[1][2];
          M0M1_2 = M[0][2]*M[1][3];
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          indx += isize_g[2];
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);

          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

          //M0M1 = M[0][2]*M[1][0];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][2]*M[1][1];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][2]*M[1][2];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];

          //M0M1 = M[0][2]*M[1][3];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];
          //indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // -//----------------------------------- //
          M0M1 = M[0][3]*M[1][0];
          M0M1_2 = M[0][3]*M[1][1];
          vM2_ = _mm256_set_ps(M[2][0], M[2][1], M[2][2], M[2][3], M[2][0], M[2][1], M[2][2], M[2][3]);
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);
          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          //----//

          M0M1 = M[0][3]*M[1][2];
          M0M1_2 = M[0][3]*M[1][3];
          vM0M1_ = _mm256_set_ps(M0M1,M0M1,M0M1,M0M1,M0M1_2,M0M1_2,M0M1_2,M0M1_2);
          vM_ = _mm256_mul_ps(vM0M1_, vM2_);

          indx += isize_g[2];
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          vPtr_ = _mm256_set_ps(ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3]);

          vVal_ = _mm256_mul_ps(vM_, vPtr_);
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

          //M0M1 = M[0][3]*M[1][0];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][3]*M[1][1];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];


          //M0M1 = M[0][3]*M[1][2];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
          //indx += isize_g[2];

          //M0M1 = M[0][3]*M[1][3];
          //ptr = &reg_grid_vals[indx + indxx];
					//val_ = M[2][0] * ptr[0];
					//val_ += M[2][1] * ptr[1];
					//val_ += M[2][2] * ptr[2];
					//val_ += M[2][3] * ptr[3];
          //val += val_ * M0M1;
		//}
          //Real val_ = sum8(vVal_);
		query_values[i] =sum8(vVal_);
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void ____vectorized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1_[2];
          Real vM[8];

          M0M1_[0] = M[0][0]*M[1][0];
          M0M1_[1] = M[0][0]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          //register Real val_;
          Real vVal[8]={0};
          Real vVal2[8]={0};
          Real* ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];


          M0M1_[0] = M[0][0]*M[1][2];
          M0M1_[1] = M[0][0]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1_[0] = M[0][1]*M[1][0];
          M0M1_[1] = M[0][1]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][1]*M[1][2];
          M0M1_[1] = M[0][1]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1_[0] = M[0][1]*M[1][0];
          M0M1_[1] = M[0][1]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][1]*M[1][2];
          M0M1_[1] = M[0][1]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];
      // ------------------------------------ //
          M0M1_[0] = M[0][3]*M[1][0];
          M0M1_[1] = M[0][3]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][3]*M[1][2];
          M0M1_[1] = M[0][3]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          val += (vVal[0]+vVal[1]+vVal[2]+vVal[3]); // * M0M1_[0];
          val += (vVal[4]+vVal[5]+vVal[6]+vVal[7]); // * M0M1_[1];
		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p



void _v2_ectorized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
  const Real lagr_denom0 = -1.0/6.0;
  const Real lagr_denom1 = 0.5;
  const Real lagr_denom2 = -0.5;
  const Real lagr_denom3 = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		point[0] = query_points[COORD_DIM * i + 0];
		grid_indx[0] = ((int)(point[0])) - 1;
		point[0] -= grid_indx[0];

		point[1] = query_points[COORD_DIM * i + 1];
		grid_indx[1] = ((int)(point[1])) - 1;
		point[1] -= grid_indx[1];

		point[2] = query_points[COORD_DIM * i + 2];
		grid_indx[2] = ((int)(point[2])) - 1;
		point[2] -= grid_indx[2];

    __m128 vx;
    vx = _mm_set1_ps(point[0]);
    __m128 c1000, c2211, c3332, vlagr;
    vlagr = _mm_set_ps(lagr_denom0,lagr_denom1,lagr_denom2,lagr_denom3);

    __m128 vM0(vlagr), vM1(vlagr), vM2(vlagr);
    __m256 vM0_0000; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM0_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM0_2222; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM0_3333; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM2_0000_1111; // elements will be M2[0] for the first 4 reg and M2[1] for the rest
    __m256 vM2_2222_3333; // elements will be M2[2] for the first 4 reg and M2[3] for the rest

    c1000 = _mm_set_ps(-1.0,-0.0,-0.0,-0.0);
    c2211 = _mm_set_ps(-2.0,-2.0,-1.0,-1.0);
    c3332 = _mm_set_ps(-3.0,-3.0,-3.0,-2.0);

    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c1000));
    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c2211));
    vM0 = _mm_mul_ps(vM0, _mm_add_ps(vx,c3332));

    vx = _mm_set1_ps(point[1]);
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c1000));
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c2211));
    vM1 = _mm_mul_ps(vM1, _mm_add_ps(vx,c3332));

    vx = _mm_set1_ps(point[2]);
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c1000));
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c2211));
    vM2 = _mm_mul_ps(vM2, _mm_add_ps(vx,c3332));
    vM2 = _mm_shuffle_ps(vM2,vM2,_MM_SHUFFLE(0, 1, 2, 3));

    vM2_0000_1111 = _mm256_set_m128(
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(3, 3, 3, 3)), // M[1][0]
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(2, 2, 2, 2)));// M[1][1]
    vM2_2222_3333 = _mm256_set_m128(
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(1, 1, 1, 1)), // M[1][2]
              _mm_shuffle_ps(vM1,vM1,_MM_SHUFFLE(0, 0, 0, 0)));// M[1][3]


    vM0_0000 = _mm256_set_m128(
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(3, 3, 3, 3)), // M[0][0]
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(3, 3, 3, 3)));// M[0][0]
    vM0_1111 = _mm256_set_m128(
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(2, 2, 2, 2)), // M[0][0]
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(2, 2, 2, 2)));// M[0][0]
    vM0_2222 = _mm256_set_m128(
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(1, 1, 1, 1)), // M[0][0]
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(1, 1, 1, 1)));// M[0][0]
    vM0_3333 = _mm256_set_m128(
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(0, 0, 0, 0)), // M[0][0]
              _mm_shuffle_ps(vM0,vM0,_MM_SHUFFLE(0, 0, 0, 0)));// M[0][0]


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
    Real* reg_ptr = &reg_grid_vals[indxx];
		Real val = 0;
    //int indx = 0;
    const int isize_g2 = isize_g[2];
    const int two_isize_g2 = 2*isize_g[2];
    const int reg_plus = isize_g[1]*isize_g2 - two_isize_g2;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          __m256 vVal_, vPtr_;
          __m256 vVal2_;
          vVal_ = _mm256_setzero_ps();


          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM2_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_0000);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM2_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_0000);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);


          reg_ptr +=  reg_plus;

      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM2_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_1111);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);


          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM2_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_1111);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);
          reg_ptr +=  reg_plus;

      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM2_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_2222);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM2_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_2222);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          reg_ptr +=  reg_plus;
      // ------------------------------------ //
          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          reg_ptr += two_isize_g2;
          vVal2_ = _mm256_mul_ps(vM2_0000_1111, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_3333);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vPtr_ = {ptr[0], ptr[1], ptr[2], ptr[3], ptr2[0], ptr2[1], ptr2[2], ptr2[3])}
          vPtr_ = _mm256_loadu2_m128(reg_ptr, reg_ptr+isize_g2);
          vVal2_ = _mm256_mul_ps(vM2_2222_3333, vPtr_);
          vVal2_ = _mm256_mul_ps(vVal2_, vM0_3333);
          vVal_ = _mm256_add_ps(vVal_, vVal2_);

          // set vm_inv = M[2][0], [1] [2] [3] in reverse order
          __m256 vM2_256 = _mm256_set_m128(vM2, vM2);
          vVal_ = _mm256_mul_ps(vVal_, vM2_256);
          val = sum8(vVal_);
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void _v1_ectorized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1 = M[0][0]*M[1][0];
          register Real val_;
          Real* ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][0]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][1]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][1]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][2]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][2]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][3]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][3]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p
#endif

void optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    Real vVal2_[8] = {0};
    Real vVal1_[8] = {0};
    int indx = 0;
    Real* ptr, *ptr2;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1;
          //register Real val_;

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][0];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][1];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][0]*vVal2_[j];


          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][2];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][3];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][0]*vVal2_[j];

          //for(int j = 0; j < 8; ++j)
          //  std::cout << "[" << j << "] = " << vVal1_[j] << std::endl;
          //do{}while(1);

          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][0];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][1];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][1]*vVal2_[j];


          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][2];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][3];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][1]*vVal2_[j];

          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][0];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][1];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][2]*vVal2_[j];

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][2];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][3];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][2]*vVal2_[j];

          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][0];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][1];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][3]*vVal2_[j];

          ptr = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          ptr2 = &reg_grid_vals[indx + indxx];
          indx += isize_g[2];
          M0M1 = M[1][2];
          vVal2_[0] = M0M1 * ptr[0];
          vVal2_[1] = M0M1 * ptr[1];
          vVal2_[2] = M0M1 * ptr[2];
          vVal2_[3] = M0M1 * ptr[3];
          M0M1 = M[1][3];
          vVal2_[4] = M0M1 * ptr2[0];
          vVal2_[5] = M0M1 * ptr2[1];
          vVal2_[6] = M0M1 * ptr2[2];
          vVal2_[7] = M0M1 * ptr2[3];
          for(int j = 0; j < 8; ++j)
            vVal1_[j] += M[0][3]*vVal2_[j];

          val = 0;
          for(int j = 0; j < 4; ++j)
            val+=vVal1_[j]*M[2][j];
          for(int j = 0; j < 4; ++j)
            val+=vVal1_[j+4]*M[2][j];
          //for(int j = 0; j < 4; ++j)
          //  std::cout << "[" << j << "] = " << vVal1_[j] * M[2][j] << std::endl;
          //for(int j = 0; j < 4; ++j)
          //  std::cout << "[" << j << "] = " << vVal1_[j+4] * M[2][j] << std::endl;

		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void gold_optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1 = M[0][0]*M[1][0];
          register Real val_;
          Real* ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][0]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][1]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][1]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][2]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][2]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][3]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][3]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void ___optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1_[2];
          Real vM[8];

          M0M1_[0] = M[0][0]*M[1][0];
          M0M1_[1] = M[0][0]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          //register Real val_;
          Real vVal[8]={0};
          Real vVal2[8]={0};
          Real* ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];


          M0M1_[0] = M[0][0]*M[1][2];
          M0M1_[1] = M[0][0]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1_[0] = M[0][1]*M[1][0];
          M0M1_[1] = M[0][1]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][1]*M[1][2];
          M0M1_[1] = M[0][1]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1_[0] = M[0][1]*M[1][0];
          M0M1_[1] = M[0][1]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][1]*M[1][2];
          M0M1_[1] = M[0][1]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];
      // ------------------------------------ //
          M0M1_[0] = M[0][3]*M[1][0];
          M0M1_[1] = M[0][3]*M[1][1];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          indx += isize_g[2];

          M0M1_[0] = M[0][3]*M[1][2];
          M0M1_[1] = M[0][3]*M[1][3];
          vM[0]=M0M1_[0];vM[1]=M0M1_[0];vM[2]=M0M1_[0];vM[3]=M0M1_[0];vM[4]=M0M1_[1];vM[5]=M0M1_[1];vM[6]=M0M1_[1];vM[7]=M0M1_[1];
          ptr = &reg_grid_vals[indx + indxx];
					vVal2[0] = M[2][0] * ptr[0];
					vVal2[1] = M[2][1] * ptr[1];
					vVal2[2] = M[2][2] * ptr[2];
					vVal2[3] = M[2][3] * ptr[3];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					vVal2[4] = M[2][0] * ptr[0];
					vVal2[5] = M[2][1] * ptr[1];
					vVal2[6] = M[2][2] * ptr[2];
					vVal2[7] = M[2][3] * ptr[3];
          for(int k = 0; k < 8; ++k)
            vVal2[k] *= vM[k];
          for(int k = 0; k < 8; ++k)
            vVal[k] += vVal2[k];
          val += (vVal[0]+vVal[1]+vVal[2]+vVal[3]); // * M0M1_[0];
          val += (vVal[4]+vVal[5]+vVal[6]+vVal[7]); // * M0M1_[1];
		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void _v4_optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
		Real val = 0;
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          Real M0M1 = M[0][0]*M[1][0];
          register Real val_;
          Real* ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][0]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][0]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][1]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][1]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][1]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][2]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][2]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][2]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];
          indx += isize_g[1]*isize_g[2] - 4 * isize_g[2];

      // ------------------------------------ //
          M0M1 = M[0][3]*M[1][0];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][1];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];


          M0M1 = M[0][3]*M[1][2];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
          indx += isize_g[2];

          M0M1 = M[0][3]*M[1][3];
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val += val_ * M0M1;
		//}
		query_values[i] = val;
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

void _optimized_interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;

	for (int i = 0; i < N_pts; i++) {
    {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
  }
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			//while (grid_indx[j] < 0)
			//	grid_indx[j] += N_reg_g[j];
		}
		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}


		const int indxx = isize_g[2] * isize_g[1] * grid_indx[0] + grid_indx[2] + isize_g[2] * grid_indx[1] ;
    register Real val_j0[4] = {0};
    int indx = 0;
		//for (int j0 = 0; j0 < 4; j0++) {
      // ------------------------------------ //
          register Real val_ = 0;
          Real* ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[0] += val_ * M[1][0];
          indx += isize_g[2];


          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[0] += val_ * M[1][1];
          indx += isize_g[2];


          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[0] += val_ * M[1][2];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[0] += val_ * M[1][3];
          indx += isize_g[1]*isize_g[2] - 3 * isize_g[2];

      // ------------------------------------ //
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[1] += val_ * M[1][0];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[1] += val_ * M[1][1];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[1] += val_ * M[1][2];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[1] += val_ * M[1][3];
          indx += isize_g[1]*isize_g[2] - 3 * isize_g[2];

      // ------------------------------------ //
          val_ = 0;
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[2] += val_ * M[1][0];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[2] += val_ * M[1][1];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[2] += val_ * M[1][2];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[2] += val_ * M[1][3];
          indx += isize_g[1]*isize_g[2] - 3 * isize_g[2];

      // ------------------------------------ //
          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[3] += val_ * M[1][0];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[3] += val_ * M[1][1];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[3] += val_ * M[1][2];
          indx += isize_g[2];

          ptr = &reg_grid_vals[indx + indxx];
					val_ = M[2][0] * ptr[0];
					val_ += M[2][1] * ptr[1];
					val_ += M[2][2] * ptr[2];
					val_ += M[2][3] * ptr[3];
          val_j0[3] += val_ * M[1][3];
          query_values[i]  = M[0][0]*val_j0[0] + M[0][1]*val_j0[1] + M[0][2]*val_j0[2] + M[0][3]*val_j0[3];
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p


/*
 * Performs a parallel 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg: The size of the original grid in each dimension.
 * @param[in] isize: The locally owned sizes that each process owns
 * @param[in] istart: The start index of each process in the global array
 * @param[in] N_pts The number of query points
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points_in The coordinates of the query points where the interpolated values are sought
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 * @param[in] c_dims: Size of the cartesian MPI communicator
 * @param[in] c_comm: MPI Communicator
 *
 */

void par_interp3_ghost_xyz_p(Real* ghost_reg_grid_vals, int data_dof,
		int* N_reg, int * isize, int* istart, const int N_pts, const int g_size,
		Real* query_points_in, Real* query_values, int* c_dims,
		MPI_Comm c_comm) {
	int nprocs, procid;
	MPI_Comm_rank(c_comm, &procid);
	MPI_Comm_size(c_comm, &nprocs);

	int N_reg_g[3], isize_g[3];
	N_reg_g[0] = N_reg[0] + 2 * g_size;
	N_reg_g[1] = N_reg[1] + 2 * g_size;
	N_reg_g[2] = N_reg[2] + 2 * g_size;

	isize_g[0] = isize[0] + 2 * g_size;
	isize_g[1] = isize[1] + 2 * g_size;
	isize_g[2] = isize[2] + 2 * g_size;

	Real h[3]; // original grid size along each axis
	h[0] = 1. / N_reg[0];
	h[1] = 1. / N_reg[1];
	h[2] = 1. / N_reg[2];

	// We copy query_points_in to query_points to aviod overwriting the input coordinates
	Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
	memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
	// Enforce periodicity
	for (int i = 0; i < N_pts; i++) {
		while (query_points[i * COORD_DIM + 0] <= -h[0]) {
			query_points[i * COORD_DIM + 0] = query_points[i * COORD_DIM + 0]
					+ 1;
		}
		while (query_points[i * COORD_DIM + 1] <= -h[1]) {
			query_points[i * COORD_DIM + 1] = query_points[i * COORD_DIM + 1]
					+ 1;
		}
		while (query_points[i * COORD_DIM + 2] <= -h[2]) {
			query_points[i * COORD_DIM + 2] = query_points[i * COORD_DIM + 2]
					+ 1;
		}

		while (query_points[i * COORD_DIM + 0] >= 1) {
			query_points[i * COORD_DIM + 0] = query_points[i * COORD_DIM + 0]
					- 1;
		}
		while (query_points[i * COORD_DIM + 1] >= 1) {
			query_points[i * COORD_DIM + 1] = query_points[i * COORD_DIM + 1]
					- 1;
		}
		while (query_points[i * COORD_DIM + 2] >= 1) {
			query_points[i * COORD_DIM + 2] = query_points[i * COORD_DIM + 2]
					- 1;
		}
	}

	// Compute the start and end coordinates that this processor owns
	Real iX0[3], iX1[3];
	for (int j = 0; j < 3; j++) {
		iX0[j] = istart[j] * h[j];
		iX1[j] = iX0[j] + (isize[j] - 1) * h[j];
	}

	// Now march through the query points and split them into nprocs parts.
	// These are stored in query_outside which is an array of vectors of size nprocs.
	// That is query_outside[i] is a vector that contains the query points that need to
	// be sent to process i. Obviously for the case of query_outside[procid], we do not
	// need to send it to any other processor, as we own the necessary information locally,
	// and interpolation can be done locally.
	int Q_local = 0, Q_outside = 0;

	// This is needed for one-to-one correspondence with output f. This is becaues we are reshuffling
	// the data according to which processor it land onto, and we need to somehow keep the original
	// index to write the interpolation data back to the right location in the output.
	std::vector<int> f_index[nprocs];
	std::vector<Real> query_outside[nprocs];
	for (int i = 0; i < N_pts; i++) {
		// The if condition check whether the query points fall into the locally owned domain or not
		if (iX0[0] - h[0] <= query_points[i * COORD_DIM + 0]
				&& query_points[i * COORD_DIM + 0] <= iX1[0] + h[0]
				&& iX0[1] - h[1] <= query_points[i * COORD_DIM + 1]
				&& query_points[i * COORD_DIM + 1] <= iX1[1] + h[1]
				&& iX0[2] - h[2] <= query_points[i * COORD_DIM + 2]
				&& query_points[i * COORD_DIM + 2] <= iX1[2] + h[2]) {
			query_outside[procid].push_back(query_points[i * COORD_DIM + 0]);
			query_outside[procid].push_back(query_points[i * COORD_DIM + 1]);
			query_outside[procid].push_back(query_points[i * COORD_DIM + 2]);
			f_index[procid].push_back(i);
			Q_local++;
			//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
			continue;
		} else {
			// If the point does not reside in the processor's domain then we have to
			// first compute which processor owns the point. After computing that
			// we add the query point to the corresponding vector.
			int dproc0 = (int) (query_points[i * COORD_DIM + 0] / h[0])
					/ isize[0];
			int dproc1 = (int) (query_points[i * COORD_DIM + 1] / h[1])
					/ isize[1];
			int proc = dproc0 * c_dims[1] + dproc1; // Compute which proc has to do the interpolation
			//PCOUT<<"proc="<<proc<<std::endl;
			query_outside[proc].push_back(query_points[i * COORD_DIM + 0]);
			query_outside[proc].push_back(query_points[i * COORD_DIM + 1]);
			query_outside[proc].push_back(query_points[i * COORD_DIM + 2]);
			f_index[proc].push_back(i);
			Q_outside++;
			//PCOUT<<"j=0 else ---------- i="<<i<<std::endl;
			continue;
		}

	}

	// Now we need to send the query_points that land onto other processor's domain.
	// This done using a sparse alltoallv.
	// Right now each process knows how much data to send to others, but does not know
	// how much data it should receive. This is a necessary information both for the MPI
	// command as well as memory allocation for received data.
	// So we first do an alltoall to get the f_index[proc].size from all processes.
	int f_index_procs_self_sizes[nprocs]; // sizes of the number of interpolations that need to be sent to procs
	int f_index_procs_others_sizes[nprocs]; // sizes of the number of interpolations that need to be received from procs

	for (int proc = 0; proc < nprocs; proc++) {
		if (!f_index[proc].empty())
			f_index_procs_self_sizes[proc] = f_index[proc].size();
		else
			f_index_procs_self_sizes[proc] = 0;
	}
	MPI_Alltoall(f_index_procs_self_sizes, 1, MPI_INT,
			f_index_procs_others_sizes, 1, MPI_INT, c_comm);

#ifdef VERBOSE2
	sleep(1);
	if(procid==0) {
		std::cout<<"procid="<<procid<<std::endl;
		std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
		std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
	}
	sleep(1);
	if(procid==1) {
		std::cout<<"procid="<<procid<<std::endl;
		std::cout<<"f_index_procs_self[0]="<<f_index_procs_self_sizes[0]<<" [1]= "<<f_index_procs_self_sizes[1]<<std::endl;
		std::cout<<"f_index_procs_others[0]="<<f_index_procs_others_sizes[0]<<" [1]= "<<f_index_procs_others_sizes[1]<<std::endl;
	}
#endif

	// Now we need to allocate memory for the receiving buffer of all query
	// points including ours. This is simply done by looping through
	// f_index_procs_others_sizes and adding up all the sizes.
	// Note that we would also need to know the offsets.
	size_t all_query_points_allocation = 0;
	int f_index_procs_others_offset[nprocs]; // offset in the all_query_points array
	int f_index_procs_self_offset[nprocs]; // offset in the query_outside array
	f_index_procs_others_offset[0] = 0;
	f_index_procs_self_offset[0] = 0;
	for (int proc = 0; proc < nprocs; ++proc) {
		// The reason we multiply by COORD_DIM is that we have three coordinates per interpolation request
		all_query_points_allocation += f_index_procs_others_sizes[proc]
				* COORD_DIM;
		if (proc > 0) {
			f_index_procs_others_offset[proc] = f_index_procs_others_offset[proc
					- 1] + f_index_procs_others_sizes[proc - 1];
			f_index_procs_self_offset[proc] =
					f_index_procs_self_offset[proc - 1]
							+ f_index_procs_self_sizes[proc - 1];
		}
	}
	int total_query_points = all_query_points_allocation / COORD_DIM;
	Real * all_query_points = (Real*) malloc(
			all_query_points_allocation * sizeof(Real));
#ifdef VERBOSE2
	if(procid==0) {
		std::cout<<"procid="<<procid<<std::endl;
		for (int proc=0;proc<nprocs;++proc)
		std::cout<<"proc= "<<proc<<" others_offset= "<<f_index_procs_others_offset[proc]<<" others_sizes= "<<f_index_procs_others_sizes[proc]<<std::endl;
		for (int proc=0;proc<nprocs;++proc)
		std::cout<<"proc= "<<proc<<" self_offset= "<<f_index_procs_self_offset[proc]<<" self_sizes= "<<f_index_procs_self_sizes[proc]<<std::endl;
	}
#endif

	MPI_Request * s_request = new MPI_Request[nprocs];
	MPI_Request * request = new MPI_Request[nprocs];

	// Now perform the allotall to send/recv query_points
	{
		int dst_r, dst_s;
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i; //(procid+i)%nprocs;
			dst_s = i; //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			int roffset = f_index_procs_others_offset[dst_r] * COORD_DIM; // notice that COORD_DIM is needed because query_points are 3 times f
			//int soffset = f_index_procs_self_offset[dst_s] * COORD_DIM;
			if (f_index_procs_others_sizes[dst_r] != 0)
				MPI_Irecv(&all_query_points[roffset],
						f_index_procs_others_sizes[dst_r] * COORD_DIM, MPI_T,
						dst_r, 0, c_comm, &request[dst_r]);
			if (!query_outside[dst_s].empty())
				MPI_Isend(&query_outside[dst_s][0],
						f_index_procs_self_sizes[dst_s] * COORD_DIM, MPI_T,
						dst_s, 0, c_comm, &s_request[dst_s]);
			//if(procid==1){
			//std::cout<<"soffset="<<soffset<<" roffset="<<roffset<<std::endl;
			//std::cout<<"f_index_procs_self_sizes[0]="<<f_index_procs_self_sizes[0]<<std::endl;
			//std::cout<<"f_index_procs_others_sizes[0]="<<f_index_procs_others_sizes[0]<<std::endl;
			//std::cout<<"q_outside["<<dst_s<<"]="<<query_outside[dst_s][0]<<std::endl;
			//}
		}
		// Wait for all the communication to finish
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
	}

	//if(procid==1){
	//  std::cout<<"total_query_points="<<total_query_points<<std::endl;
	//  std::cout<<"----- procid="<<procid<<" Q="<<all_query_points[0]<<" "<<all_query_points[1]<<" "<<all_query_points[2]<<std::endl;
	//  std::cout<<"----- procid="<<procid<<" Q="<<all_query_points[3]<<" "<<all_query_points[4]<<" "<<all_query_points[5]<<std::endl;
	//}
	//PCOUT<<"**** Q_local="<<Q_local<<" f_index_procs_self_sizes[procid]="<<f_index_procs_self_sizes[procid]<<std::endl;
	//int dum=0;
	//for (int i=0;i<Q_local;i++){
	//  dum+=query_local[i*COORD_DIM+0]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+0];
	//  dum+=query_local[i*COORD_DIM+1]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+1];
	//  dum+=query_local[i*COORD_DIM+2]-all_query_points[f_index_procs_others_offset[procid]*3+i*COORD_DIM+2];
	//}

	// Now perform the interpolation on all query points including those that need to
	// be sent to other processors and store them into all_f_cubic
	Real* all_f_cubic = (Real*) malloc(
			total_query_points * sizeof(Real) * data_dof);
	interp3_ghost_xyz_p(ghost_reg_grid_vals, data_dof, N_reg, N_reg_g, isize_g,
			istart, total_query_points, g_size, all_query_points, all_f_cubic);

	//if(procid==0){
	//  std::cout<<"total_query_points="<<total_query_points<<std::endl;
	//  std::cout<<"procid="<<procid<<" Q="<<all_query_points[0]<<" "<<all_query_points[1]<<" "<<all_query_points[2]<<" f= "<<all_f_cubic[0]<<std::endl;
	//  std::cout<<"procid="<<procid<<" Q="<<all_query_points[3]<<" "<<all_query_points[4]<<" "<<all_query_points[5]<<" f= "<<all_f_cubic[1]<<std::endl;
	//}

	// Now we have to do an alltoall to distribute the interpolated data from all_f_cubic to
	// f_cubic_unordered.
	Real * f_cubic_unordered = (Real*) malloc(N_pts * sizeof(Real) * data_dof); // The reshuffled semi-final interpolated values are stored here
	{
		//PCOUT<<"total_query_points="<<total_query_points<<" N_pts="<<N_pts<<std::endl;
		int dst_r, dst_s;
		MPI_Datatype stype[nprocs], rtype[nprocs];
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_vector(data_dof, f_index_procs_self_sizes[i], N_pts, MPI_T,
					&rtype[i]);
			MPI_Type_vector(data_dof, f_index_procs_others_sizes[i],
					total_query_points, MPI_T, &stype[i]);
			MPI_Type_commit(&stype[i]);
			MPI_Type_commit(&rtype[i]);
		}
		for (int i = 0; i < nprocs; ++i) {
			dst_r = i; //(procid+i)%nprocs;
			dst_s = i; //(procid-i+nprocs)%nprocs;
			s_request[dst_s] = MPI_REQUEST_NULL;
			request[dst_r] = MPI_REQUEST_NULL;
			// Notice that this is the adjoint of the first comm part
			// because now you are sending others f and receiving your part of f
			int soffset = f_index_procs_others_offset[dst_r];
			int roffset = f_index_procs_self_offset[dst_s];
			//if(procid==0)
			//  std::cout<<"procid="<<procid<<" dst_s= "<<dst_s<<" soffset= "<<soffset<<" s_size="<<f_index_procs_others_sizes[dst_s]<<" dst_r= "<<dst_r<<" roffset="<<roffset<<" r_size="<<f_index_procs_self_sizes[dst_r]<<std::endl;
			//if(f_index_procs_self_sizes[dst_r]!=0)
			//  MPI_Irecv(&f_cubic_unordered[roffset],f_index_procs_self_sizes[dst_r],rtype, dst_r,
			//      0, c_comm, &request[dst_r]);
			//if(f_index_procs_others_sizes[dst_s]!=0)
			//  MPI_Isend(&all_f_cubic[soffset],f_index_procs_others_sizes[dst_s],stype,dst_s,
			//      0, c_comm, &s_request[dst_s]);
			//
			if (f_index_procs_self_sizes[dst_r] != 0)
				MPI_Irecv(&f_cubic_unordered[roffset], 1, rtype[i], dst_r, 0,
						c_comm, &request[dst_r]);
			if (f_index_procs_others_sizes[dst_s] != 0)
				MPI_Isend(&all_f_cubic[soffset], 1, stype[i], dst_s, 0, c_comm,
						&s_request[dst_s]);
		}
		MPI_Status ierr;
		for (int proc = 0; proc < nprocs; ++proc) {
			if (request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&request[proc], &ierr);
			if (s_request[proc] != MPI_REQUEST_NULL)
				MPI_Wait(&s_request[proc], &ierr);
		}
		for (int i = 0; i < nprocs; ++i) {
			MPI_Type_free(&stype[i]);
			MPI_Type_free(&rtype[i]);
		}
	}

	// Now copy back f_cubic_unordered to f_cubic in the correct f_index
	for (int dof = 0; dof < data_dof; ++dof) {
		for (int proc = 0; proc < nprocs; ++proc) {
			if (!f_index[proc].empty())
				for (int i = 0; i < (int)f_index[proc].size(); ++i) {
					int ind = f_index[proc][i];
					//f_cubic[ind]=all_f_cubic[f_index_procs_others_offset[proc]+i];
					query_values[ind + dof * N_pts] =
							f_cubic_unordered[f_index_procs_self_offset[proc]
									+ i + dof * N_pts];
				}
		}
	}

	free(query_points);
	free(all_query_points);
	free(all_f_cubic);
	free(f_cubic_unordered);
	delete[] s_request;
	delete[] request;
	//vector
	for (int proc = 0; proc < nprocs; ++proc) {
		std::vector<int>().swap(f_index[proc]);
		std::vector<Real>().swap(query_outside[proc]);
	}
	return;
} // end of par_interp3_ghost_xyz_p

// the factor is computed by the following transform:
// X=0 -> Y = ghp
// X=1-h -> Y = 1-hp - ghp
void rescale(const int g_size, int* N_reg, int* N_reg_g, int* istart,
		const int N_pts, Real* query_points) {

	if (g_size == 0)
		return;
	Real hp[3];
	Real h[3];
	hp[0] = 1. / N_reg_g[0]; // New mesh size
	hp[1] = 1. / N_reg_g[1]; // New mesh size
	hp[2] = 1. / N_reg_g[2]; // New mesh size

	h[0] = 1. / (N_reg[0]); // old mesh size
	h[1] = 1. / (N_reg[1]); // old mesh size
	h[2] = 1. / (N_reg[2]); // old mesh size

	Real factor[3];
	factor[0] = (1. - (2. * g_size + 1.) * hp[0]) / (1. - h[0]);
	factor[1] = (1. - (2. * g_size + 1.) * hp[1]) / (1. - h[1]);
	factor[2] = (1. - (2. * g_size + 1.) * hp[2]) / (1. - h[2]);
	for (int i = 0; i < N_pts; i++) {
		query_points[0 + COORD_DIM * i] = (query_points[0 + COORD_DIM * i]
				- istart[0] * h[0]) * factor[0] + g_size * hp[0];
		query_points[1 + COORD_DIM * i] = (query_points[1 + COORD_DIM * i]
				- istart[1] * h[1]) * factor[1] + g_size * hp[1];
		//query_points[2+COORD_DIM*i]=(query_points[2+COORD_DIM*i]-istart[2]*h[2])*factor[2]+g_size*hp[2];
	}
	return;
} // end of rescale

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 * snafu
 *
 */

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[4];
  lagr_denom[0] = -1.0/6.0;
  lagr_denom[1] = 0.5;
  lagr_denom[2] = -0.5;
  lagr_denom[3] = 1.0/6.0;
	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j];// * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						//val += M[0][j0] * M[1][j1] * M[2][j2]
						//		* reg_grid_vals[0];
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			//query_values[0] = val;
			query_values[i + k * N_pts] = val;
		}
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on all sides
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[in] interp_order The order of interpolation (e.g. 3 for cubic)
 * @param[out] query_values The interpolated values
 *
 */

void interp3_ghost_xyz_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values,
    int interp_order,
		bool query_values_already_scaled) {
	Real* query_points;

	if (query_values_already_scaled == false) {
		// First we need to rescale the query points to the new padded dimensions
		// To avoid changing the user's input we first copy the query points to a
		// new array
		query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
		memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
		rescale_xyz(g_size, N_reg, N_reg_g, istart, N_pts, query_points);
	} else {
		query_points = query_points_in;
	}
	Real lagr_denom[interp_order + 1];
	for (int i = 0; i < interp_order + 1; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < interp_order + 1; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];

	for (int i = 0; i < N_pts; i++) {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];

		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}

#ifdef VERBOSE2
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
#endif

		Real M[3][interp_order + 1];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < interp_order + 1; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < interp_order + 1; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < interp_order + 1; j2++) {
				for (int j1 = 0; j1 < interp_order + 1; j1++) {
					for (int j0 = 0; j0 < interp_order + 1; j0++) {
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}

	if (query_values_already_scaled == false) {
		free(query_points);
	}
	return;

}  // end of interp3_ghost_xyz_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * This function assumes that the input grid values have been padded on x, and y directions
 * by g_size grids.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] g_size The number of ghost points padded around the input array
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_ghost_p(Real* reg_grid_vals, int data_dof, int* N_reg,
		int* N_reg_g, int * isize_g, int* istart, const int N_pts,
		const int g_size, Real* query_points_in, Real* query_values) {

	// First we need to rescale the query points to the new padded dimensions
	// To avoid changing the user's input we first copy the query points to a
	// new array
	Real* query_points = (Real*) malloc(N_pts * COORD_DIM * sizeof(Real));
	memcpy(query_points, query_points_in, N_pts * COORD_DIM * sizeof(Real));
	rescale(g_size, N_reg, N_reg_g, istart, N_pts, query_points);

	//std::cout<<"N_reg[0]="<<N_reg[0]<<" N_reg[1]="<<N_reg[1]<<" N_reg[2]="<<N_reg[2]<<std::endl;
	//std::cout<<"N_reg_g[0]="<<N_reg_g[0]<<" N_reg_g[1]="<<N_reg_g[1]<<" N_reg_g[2]="<<N_reg_g[2]<<std::endl;

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = isize_g[0] * isize_g[1] * isize_g[2];
	//int N_pts=query_points.size()/COORD_DIM;

	for (int i = 0; i < N_pts; i++) {
#ifdef VERBOSE2
		std::cout<<"q[0]="<<query_points[i*3+0]<<std::endl;
		std::cout<<"q[1]="<<query_points[i*3+1]<<std::endl;
		std::cout<<"q[2]="<<query_points[i*3+2]<<std::endl;
#endif
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		//grid_indx[0]=15;
		//grid_indx[1]=15;
		//grid_indx[2]=14;
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg_g[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg_g[j];
		}
#ifdef VERBOSE2
		std::cout<<"***** grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;
		std::cout<<"***** point="<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
		std::cout<<"f @grid_index="<<reg_grid_vals[grid_indx[0]*isize_g[1]*isize_g[2]+grid_indx[1]*isize_g[2]+grid_indx[2]]<<std::endl;
		std::cout<<"hp= "<<1./N_reg_g[0]<<std::endl;
		std::cout<<"N_reg_g= "<<N_reg_g[0]<<" "<<N_reg_g[1]<<" "<<N_reg_g[2]<<std::endl;
#endif

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[2]+j2)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[0]+j0)%N_reg);
						//int indx = ((grid_indx[2]+j2)%N_reg_g[2]) + N_reg_g[2]*((grid_indx[1]+j1)%N_reg_g[1]) + N_reg_g[2]*N_reg_g[1]*((grid_indx[0]+j0)%N_reg_g[0]);
						int indx = ((grid_indx[2] + j2) % isize_g[2])
								+ isize_g[2]
										* ((grid_indx[1] + j1) % isize_g[1])
								+ isize_g[2] * isize_g[1]
										* ((grid_indx[0] + j0) % isize_g[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
	free(query_points);
}            //end of interp3_ghost_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) )
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg An integer pointer that specifies the size of the grid in each dimension.
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_p(Real* reg_grid_vals, int data_dof, int* N_reg, const int N_pts,
		Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg[0] * N_reg[1] * N_reg[2];
	//int N_pts=query_points.size()/COORD_DIM;
	//query_values.resize(N_pts*data_dof);

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg[j];
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg[j];
		}
		//std::cout<<"grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[2]+j2)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[0]+j0)%N_reg);
						int indx = ((grid_indx[2] + j2) % N_reg[2])
								+ N_reg[2] * ((grid_indx[1] + j1) % N_reg[1])
								+ N_reg[2] * N_reg[1]
										* ((grid_indx[0] + j0) % N_reg[0]);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p

/*
 * Performs a 3D cubic interpolation for a row major periodic input (x \in [0,1) ).
 * Limitation: The grid must be cubic, i.e. the number of grid points must be the same
 * in all dimensions.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg The size of the regular grid (The grid must have the same size in all dimensions)
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */

void interp3_p(Real* reg_grid_vals, int data_dof, int N_reg, const int N_pts,
		Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg * N_reg * N_reg;
	//int N_pts=query_points.size()/COORD_DIM;
	//query_values.resize(N_pts*data_dof);

	for (int i = 0; i < N_pts; i++) {
		Real point[COORD_DIM];
		int grid_indx[COORD_DIM];
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg;
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg;
		}
		//std::cout<<"grid_index="<<grid_indx[0]<<" "<<grid_indx[1]<<" "<<grid_indx[2]<<std::endl;

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						//int indx = ((grid_indx[0]+j0)%N_reg) + N_reg*((grid_indx[1]+j1)%N_reg) + N_reg*N_reg*((grid_indx[2]+j2)%N_reg);
						int indx = ((grid_indx[2] + j2) % N_reg)
								+ N_reg * ((grid_indx[1] + j1) % N_reg)
								+ N_reg * N_reg * ((grid_indx[0] + j0) % N_reg);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p

/*
 * Performs a 3D cubic interpolation for a column major periodic input (x \in [0,1) )
 * Limitation: The grid must be cubic, i.e. the number of grid points must be the same
 * in all dimensions.
 * @param[in] reg_grid_vals The function value at the regular grid
 * @param[in] data_dof The degrees of freedom of the input function. In general
 * you can input a vector to be interpolated. In that case, each dimension of the
 * vector should be stored in a linearized contiguous order. That is the first dimension
 * should be stored in reg_grid_vals and then the second dimension, ...
 *
 * @param[in] N_reg The size of the regular grid (The grid must have the same size in all dimensions)
 *
 * @param[in] N_pts The number of query points
 *
 * @param[in] query_points The coordinates of the query points where the interpolated values are sought.
 * One must store the coordinate values back to back. That is each 3 consecutive values in its array
 * determine the x,y, and z coordinate of 1 query point in 3D.
 *
 * @param[out] query_values The interpolated values
 *
 */
void interp3_p_col(Real* reg_grid_vals, int data_dof, int N_reg,
		const int N_pts, Real* query_points, Real* query_values) {

	Real lagr_denom[4];
	for (int i = 0; i < 4; i++) {
		lagr_denom[i] = 1;
		for (int j = 0; j < 4; j++) {
			if (i != j)
				lagr_denom[i] /= (Real) (i - j);
		}
	}

	int N_reg3 = N_reg * N_reg * N_reg;

	Real point[COORD_DIM];
	int grid_indx[COORD_DIM];
	for (int i = 0; i < N_pts; i++) {
		for (int j = 0; j < COORD_DIM; j++) {
			point[j] = query_points[COORD_DIM * i + j] * N_reg;
			grid_indx[j] = (floor(point[j])) - 1;
			point[j] -= grid_indx[j];
			while (grid_indx[j] < 0)
				grid_indx[j] += N_reg;
		}

		Real M[3][4];
		for (int j = 0; j < COORD_DIM; j++) {
			Real x = point[j];
			for (int k = 0; k < 4; k++) {
				M[j][k] = lagr_denom[k];
				for (int l = 0; l < 4; l++) {
					if (k != l)
						M[j][k] *= (x - l);
				}
			}
		}

		for (int k = 0; k < data_dof; k++) {
			Real val = 0;
			for (int j2 = 0; j2 < 4; j2++) {
				for (int j1 = 0; j1 < 4; j1++) {
					for (int j0 = 0; j0 < 4; j0++) {
						int indx = ((grid_indx[0] + j0) % N_reg)
								+ N_reg * ((grid_indx[1] + j1) % N_reg)
								+ N_reg * N_reg * ((grid_indx[2] + j2) % N_reg);
						val += M[0][j0] * M[1][j1] * M[2][j2]
								* reg_grid_vals[indx + k * N_reg3];
					}
				}
			}
			query_values[i + k * N_pts] = val;
		}
	}
} // end of interp3_p_col

