#pragma once
#include <cassert>
#include <omp.h>

namespace converter
{
  namespace dim_2
  {
    //Define some Physical constants
    constexpr double C_LIGHT=299792458.; // [m/s] 

    
    extern "C" {
      struct part{
	int   l_x, l_y, l_w, l_px, l_py, l_pz, l_vx, l_vy, l_vz, l_ek, l_rm, l_gm;
	double *x, *y, *w,  *px, *py, *pz, *vx, *vy, *vz, *ek, *rm, *gm;
      };
      
      void read_particle(const char* cstring, int clen,
			 const char* spec,    int len_sp,
			 part* arrays, long* npart, long* npart_proc, long* start);
      
      struct field{
	int   global_sx, global_sy, l_sx, l_sy, l_dx, l_dy, l_gridx, l_gridy, l_data_size;
	double stagger;
	double* gridx,  *gridy, *l_field_data;
      };
      
      void read_fields( const char* cstring, int clen,  const char* blockid, int blen, field* field_x);
      void read_derived_vars( const char* cstring, int clen,  const char* blockid, int blen, field* field_x);        
      void init_read();
      
    }
       
    static inline bool sort_using_less_than(double u, double v)
    {
      return u < v;
    }
    
    static inline double  max(double A[], int len)
    {
      int max_val = A[0];
 
      #pragma omp parallel for reduction(max:max_val) 
      for (int idx = 0; idx < len; idx++)
	max_val = max_val > A[idx] ? max_val : A[idx];
      
      return max_val;
    }

    static inline double  min(double A[], int len)
    {
      int min_val = A[0];  

      #pragma omp parallel for reduction(min:min_val) 
      for (int idx = 0; idx < len; idx++)
	min_val = min_val < A[idx] ? min_val : A[idx];
      
      return min_val;
    }



    template<typename T> 
    static bool isLess(T left, T right, bool orequal = false) {
      T eps = static_cast<T>(std::numeric_limits<T>::epsilon);
      if (std::fabs(left - right) < eps){
	return (orequal);
      }
      return (left < right);
    }
  
    template<typename T> 
    static bool isGreater(T left, T right, bool orequal = false) {
      if (std::fabs(left - right) < std::numeric_limits<T>::epsilon) {
	return (orequal);
      }
      return (left > right);
    }
    
    
    // Value is within the range [low, high).
    template <typename T>      
      bool is_inside(const T& value, const T& low, const T& high) {
      if ((high==0) && (low==0)) return true;
      if (high==low) return true;
      return !(value < low) && (value < high);
      //return isGreater<T>(value,low) && isLess<T>(value, high);
    }


    double re_binning(const std::vector<double> &p){
      // Re-binning using Freedman-Diaconis' rule.
      size_t n = p.size();
      double prange = *std::max_element(p.begin(), p.end()) -
	*std::min_element(p.begin(), p.end());
      double bin_w = 1.0;
      if (n>1) {
	size_t q1  = static_cast<size_t>(n * 0.25);
	size_t q3 = n - static_cast<size_t>(n * 0.25);
	auto p_c = p;
	std::nth_element(p_c.begin(), p_c.begin() + q1, p_c.end());
	std::nth_element(p_c.begin(), p_c.begin() + q3, p_c.end());
	double iq_r = p_c[q3] - p_c[q1];
	double iq = std::max(iq_r, prange / 10.);
	bin_w = 2 * iq * pow(n, -1. / 3.);
      }
      return bin_w;
    }
  
    
    
  }//!dim_2
}//!converter
