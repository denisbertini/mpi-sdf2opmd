#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>

// MPI
#include "mpi.h"

// Utils
#include "utils.h"

using namespace std;

namespace converter
{
  namespace dim_3
  {
    
    class Merger_3d {
    private:
      int m_verbosity{0};
      int m_min_npart_pcell{0};
      double m_min_dp_pcell[3];
      int m_mpi_rank{-1};
      int m_mpi_size{-1};
      std::vector<int> m_mask_indexes;
      bool m_auto_binning{false};
      int m_total_part_processed{0};      
    public:
      Merger_3d(int rank, int size);
      virtual ~Merger_3d(){;}
      void merge(part &pp, const int* x_bins, const int* p_bins);
      void p_cartesian(part &pp, const int* x_bins, const int* p_bins, std::vector<int>* c_indexes);
      std::tuple<int, int> p_reduction(part &pp, const int* p_bins, const double* p_min,
				       const double* dp, std::vector<int>* c_indexes);
      void setVerbose(int val){m_verbosity=val;}
      void setMinNpartPerCell(int val){m_min_npart_pcell=val;}
      void setMinDpPerCell(int i, double val){m_min_dp_pcell[i]=val;}
      void setMpiInfo(int i, int j){m_mpi_rank=i;m_mpi_size=j;}
      std::vector<int> get_mask_indexes(){return m_mask_indexes;}
    };
  }
}
