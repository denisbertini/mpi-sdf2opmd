#include "merger_3d.h"
#include <chrono>
#include <omp.h>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

namespace converter
{
  namespace dim_3
  {
  
    Merger_3d::Merger_3d(int rank, int size)
      : m_verbosity(0)			   
      , m_min_npart_pcell(4)
      , m_min_dp_pcell{1e-10, 1e-10, 1e-10} 
    {
      m_mpi_rank=rank;
      m_mpi_size=size;
    }
    
    void Merger_3d::merge(part &pp, const int* x_bins, const int* p_bins){    

      //  
      // Particle resampling, Vranic and al, method
      //
      
      if (m_verbosity>1 && 0 == m_mpi_rank){
	std::cout << std::endl;
	std::cout << "Merger: rank: "<< m_mpi_rank
		  << " (x,y,z)   binning -> (" << x_bins[0] << "," << x_bins[1] << "," << x_bins[2] <<  ")"
		  << " (px,py,pz) binning -> (" << p_bins[0] << "," << p_bins[1] << "," << p_bins[2] << ")"
	          << " part: len_z: " << pp.l_z 
		  << std:: endl;
	std::cout << std::endl;
      }

       // Get parallel min max
      auto t1 = high_resolution_clock::now();
      
      double xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0;
      double zmin=0.0, zmax=0.0;
      xmax = max(pp.x, pp.l_x);
      ymax = max(pp.y, pp.l_y);
      zmax = max(pp.z, pp.l_z);
      xmin = min(pp.x, pp.l_x); 
      ymin = min(pp.y, pp.l_y);
      zmin = min(pp.z, pp.l_z);      
      
      auto t2 = high_resolution_clock::now();
      
      duration<double, std::milli> ms_double1 = t2 - t1;
      if ((m_verbosity>1) && (0 == m_mpi_rank)){ 
	std::cout << " Find min-max timing: "  << ms_double1.count() << "ms\n"
		  << std::endl;     	
	std::cout << " xmin: " << xmin << " ymin: " << ymin << " zmin: " <<  zmin << std::endl; 
	std::cout << " xmax: " << xmax<< "  ymax: " << ymax << " zmax: " <<  zmax << std::endl; 
      }

      int npart=pp.l_x;
      int dims = x_bins[0]*x_bins[1]*x_bins[2]; 

      double delta_x = fabs(xmax-xmin)/x_bins[0];
      double delta_y = fabs(ymax-ymin)/x_bins[1];
      double delta_z = fabs(zmax-zmin)/x_bins[2];
      
      std::vector<int> cells_indexes[dims];
             
      for (int i=0;i<x_bins[0];i++){
	for(int j=0;j<x_bins[1]; j++){
	  for(int k=0; k<x_bins[2]; k++){
	    int  cell_index  = ( i * x_bins[1] + j ) * x_bins[2] + k;

	      if ((m_verbosity>1) && (0 == m_mpi_rank)){ 
		std::cout << "3D linear cell index i: " << i
			  << " j: " << j  << " k: " << k 
			  << " lin: " << cell_index
			  << std::endl;
		int i_x = (int) ((cell_index)/(x_bins[2]*x_bins[1])); 
		int i_y = (int) (cell_index/x_bins[2]) % x_bins[1];
		int i_z = (int) (cell_index % x_bins[2]);
		std::cout << "3D recalcuted  cell_index--> "
			  <<  " i_x: " << i_x 
			  <<  " i_y: " << i_y
			  <<  " i_z: " << i_z
			  << std::endl;
	      }
	      
	      double x_inf = xmin + i * delta_x;
	      double x_sup = xmin + (i+1) * delta_x;
	      double y_inf = ymin + j * delta_y;
	      double y_sup = ymin + (j+1) * delta_y;
	      double z_inf = zmin + k* delta_z;
	      double z_sup = zmin + (k+1) * delta_z;	      
	      
	      for (int l=0; l<npart;l++){
		if ( is_inside<double>(pp.x[l], x_inf, x_sup) &&
		     is_inside<double>(pp.y[l], y_inf, y_sup) &&
		     is_inside<double>(pp.z[l], z_inf, z_sup) 
		     )
		cells_indexes[cell_index].push_back(l);
	      }//!l
	    }//!k
	  }//!j
	}//!i

	

      int sum_points = 0;      
      if ((m_verbosity>1) && (0==m_mpi_rank)){
	for(int  i = 0; i < dims ; ++i){
	  int i_x = (int) ((i)/(x_bins[2]*x_bins[1])); 
	  int i_y = (int) (i/x_bins[2]) % x_bins[1];
	  int i_z = (int) (i % x_bins[2]);

	  std::cout << " cell_index : "<< i
		    << " x[" << xmin +i_x*delta_x  << "," << xmin+(i_x+1)*delta_x << "["
		    << " y[" << ymin+ i_y*delta_y << "," << ymin+ (i_y+1)*delta_y << "["
		    << " z[" << zmin+ i_z*delta_z << "," << zmin+ (i_z+1)*delta_z << "["	    	    
	            <<" contains npart:  "
		    << cells_indexes[i].size() << std::endl;
	  
	  sum_points+=cells_indexes[i].size(); 	   
	}
      }

      // Cross check mapping
      assert(sum_points != pp.l_x && "non-consistent particle to mesh (x) mapping!");
      assert(sum_points != pp.l_y && "non-consistent particle to mesh (y) mapping!");	
      assert(sum_points != pp.l_z && "non-consistent particle to mesh (z) mapping!");	

      // Cartesian momentum merging using Vranic & al. method  
      p_cartesian(pp, x_bins, p_bins, cells_indexes);
      
    }//! merge


    
    void  Merger_3d::p_cartesian(part &pp, const int* x_bins, const int* p_bins,
				 std::vector<int> c_indexes[]){
      //
      // Momentum cell analysis
      //

      if ((m_verbosity>1) && (0 == m_mpi_rank)){  
      std::cout << "p_cartesian rank: " << m_mpi_rank << std::endl;
	for(int  i = 0; i < x_bins[0]*x_bins[1]*x_bins[2] ; ++i){
	  std::cout << " cell_index: "
		    << i <<" contains npart: "
		    << c_indexes[i].size() << std::endl;
	}		
      }

      int dims= x_bins[0]*x_bins[1]*x_bins[2];
      int pdims=p_bins[0]*p_bins[1]*p_bins[2];

      //Main loop over linearized cell indexes
      for(int  c_index = 0; c_index < dims ; ++c_index){
	// Check nb. of particles in the cell
	int npart_pcell = c_indexes[c_index].size();	
	if (npart_pcell < m_min_npart_pcell) continue;	

	// Define momentum tiles per cell 
	std::vector<double> p_cell[3];
	for(int i=0;i<3;i++) p_cell[i].resize(npart_pcell);
        for(int ic=0;ic<npart_pcell; ic++){
	  p_cell[0][ic] = pp.px[c_indexes[c_index][ic]];
	  p_cell[1][ic] = pp.py[c_indexes[c_index][ic]];
	  p_cell[2][ic] = pp.pz[c_indexes[c_index][ic]];    	    
        }

	
	if ((m_verbosity>1) && (0 == m_mpi_rank)){
	  std::cout << " icell Linear: " << c_index
	       	    << " npart/cell:" << npart_pcell
		    << std::endl;
	  
	   for(int ii=0;ii<npart_pcell;ii++) std::cout 
					      << " px: " << p_cell[0][ii] 
					      << " py: " << p_cell[1][ii]
					      << " pz: " << p_cell[2][ii]
					      << std::endl;					      	  	  
	}
	
	//Define new binning in  momentum space  
	double p_min[3], p_max[3], dp[3], inv_dp[3];
	bool p_binning[3]={true,true,true};
	double p_rebin[3]={0.,0.,0.}; 
	int n_rebin[3]={0,0,0};
	
	for(int i=0;i<3;i++){
	  auto mm = minmax_element(p_cell[i].begin(), p_cell[i].end());
	  p_min[i] = *mm.first;
	  p_max[i] = *mm.second;
	  // Mark anomaly in direction space
	  if (fabs(p_max[i]-p_min[i]) == 0) p_binning[i] = false;	  
	  p_rebin[i] = re_binning(p_cell[i]);
	  //n_rebin[i] = (p_rebin[i]==0) ? 0 : std::ceil((p_max[i]-p_min[i])/p_rebin[i]);
	  n_rebin[i] = (p_rebin[i]==0) ? 1 : std::ceil((p_max[i]-p_min[i])/p_rebin[i]);
	  if (p_binning[i]){
	    dp[i] = fabs(p_max[i] - p_min[i])/(n_rebin[i]);
	    inv_dp[i] = 1.0/dp[i];
	  }else{
	    dp[i]=0.0;
	    inv_dp[i]=1.0;
	  }
	  
	  if ((m_verbosity>1) && (0 == m_mpi_rank)){
	    std::cout << " icell Linear: " << c_index
	              << " dir= " << i
		      << " pmin: " << p_min[i]
		      << " pmax: " << p_max[i]
		      << " n_pbins: " << n_rebin[i]
	              << " inv_dp: " << inv_dp[i] 
	              << " p_binning: " << p_binning[i] 
		      << std::endl;	    
	  }	  
	}//!for(3_D)

	// Do the mapping p_cell <-> epoch_trk_indexes
	pdims=n_rebin[2]*n_rebin[1]*n_rebin[0];
	std::vector<int> p_cells_indexes[pdims];	
	for (int i=0;i<n_rebin[0];i++){
	  for(int j=0;j<n_rebin[1]; j++){
	    for(int k=0;k<n_rebin[2]; k++){
	      //Row major mapping 
	      int p_cell_index = (i*n_rebin[1]*n_rebin[2])+(j*n_rebin[2])+k;
	      
	      if ((m_verbosity>1) && (0 == m_mpi_rank)){ 
		std::cout << "3D linear cell index i: " << i
			  << " j: " << j  << " k: " << k 
			  << " lin: " << p_cell_index
			  << std::endl;
		int i_px = (int) ((p_cell_index)/(n_rebin[2]*n_rebin[1])); 
		int i_py = (int) (p_cell_index/n_rebin[2]) % n_rebin[1];
		int i_pz = (int) (p_cell_index % n_rebin[2]);
		std::cout << "3D recalcuted  p_index--> "
			  <<  " i_px: " << i_px 
			  <<  " i_py: " << i_py
			  <<  " i_pz: " << i_pz
			  << std::endl;
	      }
	      
	      double px_inf = p_min[0] + i*dp[0];
	      double px_sup = p_min[0] + (i+1) * dp[0];
	      double py_inf = p_min[1] + j* dp[1];
	      double py_sup = p_min[1] + (j+1) * dp[1];
	      double pz_inf = p_min[2] + k* dp[2];
	      double pz_sup = p_min[2] + (k+1) * dp[2];	      
	      
	      for (int l=0; l<npart_pcell;l++){
		if ( is_inside<double>(p_cell[0][l], px_inf, px_sup) &&
		     is_inside<double>(p_cell[1][l], py_inf, py_sup) &&
		     is_inside<double>(p_cell[2][l], pz_inf, pz_sup) 
		     )
		  p_cells_indexes[p_cell_index].push_back(c_indexes[c_index][l]);
	      }//!l
	    }//!k
	  }//!j
	}//!i
	
	// Do particle reduction
	p_reduction(pp,n_rebin,p_min,dp,p_cells_indexes);
	
	// Check particle distribution after Mapping 
	for (int i_pcell=0; i_pcell<pdims; i_pcell++){
	  if ((m_verbosity>1) && (0 == m_mpi_rank)){
	    if (p_cells_indexes[i_pcell].size()==0) continue;
	    std::cout << "Particle distribution: i_pcell: " << i_pcell << " npart/cell: "
		      << p_cells_indexes[i_pcell].size() << std::endl;
	    for (size_t p=0; p<p_cells_indexes[i_pcell].size(); p++)
	      std::cout << " p_index: " << p
			<< " px: " << pp.px[p_cells_indexes[i_pcell][p]]
			<< " py: " << pp.py[p_cells_indexes[i_pcell][p]]
			<< " pz: " << pp.pz[p_cells_indexes[i_pcell][p]]
			<< std::endl;
	  }
	}
	
      }//!for(c_index)		      
    }//!p_cartesian

    
    void Merger_3d::p_reduction(part &pp, const int* p_bins, const double* p_min,
				const double* dp, std::vector<int> p_indexes[]){

      //
      // Particle Merging Algorithm
      //          M. Vranic et al., CPC, 191 65-73 (2015)
      //          inspired by specific SMILEI implementation
      //          https://smileipic.github.io/Smilei/Understand/particle_merging.html 
      //
      
      // Photon case not treated for the moment 
      if (!p_indexes[0].empty() && pp.rm[p_indexes[0][0]]==0) return;
      
      // P-Cell geometry
      int pdims{0};       
      // Total number of P-bins
      pdims = p_bins[0]*p_bins[1]*p_bins[2];

      // Check particle distribution after Mapping 
      for (int i_pcell=0; i_pcell<pdims; i_pcell++){
	int npart_per_cell=p_indexes[i_pcell].size();
	
	// Selected only enough populated cells.
	if (npart_per_cell<m_min_npart_pcell) continue;
	
	if ((m_verbosity>1) && (0 == m_mpi_rank)){	 
	  std::cout << "p_reduction: selected i_pcell: " << i_pcell << " npart/cell: "
		    << npart_per_cell << std::endl;
	}
	
	// Init total cell quantities
	double tot_w{0.0};
	double tot_px{0.0};
	double tot_py{0.0};
	double tot_pz{0.0};
	double tot_e{0.0};   
	double mo{0.0};
	
	for (int p=0; p<npart_per_cell; p++)
	  {
	    // Total weight (wt)
	    tot_w  += pp.w[p_indexes[i_pcell][p]];	    
	    // total momentum  (pt)
	    tot_px += pp.px[p_indexes[i_pcell][p]]*pp.w[p_indexes[i_pcell][p]];
	    tot_py += pp.py[p_indexes[i_pcell][p]]*pp.w[p_indexes[i_pcell][p]];
	    tot_pz += pp.pz[p_indexes[i_pcell][p]]*pp.w[p_indexes[i_pcell][p]];	    
	    // total energy
	    tot_e +=  pp.rm[p_indexes[i_pcell][p]]*pp.w[p_indexes[i_pcell][p]]*C_LIGHT*C_LIGHT;
	    mo = (mo==0)? pp.rm[p_indexes[i_pcell][p]]/pp.gm[p_indexes[i_pcell][p]] : mo;
	    
	    if ((m_verbosity>1) && (0 == m_mpi_rank)){	 
	      std::cout.precision(8);
	      std::cout << " p_index: " << p
			<< " px: " << pp.px[p_indexes[i_pcell][p]] << " [kg m s^-1] "
			<< " py: " << pp.py[p_indexes[i_pcell][p]] << " [kg m s^-1] "
			<< " pz: " << pp.pz[p_indexes[i_pcell][p]] << " [kg m s^-1] "		  
			<< " gm: " << pp.gm[p_indexes[i_pcell][p]] << " [kg m s^-1] "
			<< " rm: " << pp.rm[p_indexes[i_pcell][p]] << " [kg] "
			<< " me: " << pp.rm[p_indexes[i_pcell][p]]/pp.gm[p_indexes[i_pcell][p]] << " [kg] "
			<< " Ee: " << pp.rm[p_indexes[i_pcell][p]]*C_LIGHT*C_LIGHT << " [J] (kg m^2 c^-2) "
			<< std::endl;
	    }
	  }//!for(particle/cell)
	
	// 3D index map <check-me>
	int i_px = (int) ((i_pcell)/(p_bins[2]*p_bins[1])); 
	int i_py = (int) (i_pcell/p_bins[2]) % p_bins[1];
	int i_pz = (int) (i_pcell % p_bins[2]);
	
	if ((m_verbosity>1) && (0 == m_mpi_rank)){	 
	  std::cout << "P-cell dim: " << p_bins[0]*p_bins[1]*p_bins[2]
		    << " linear: " << i_pcell 
		    << " i_px: " << i_px
		    << " i_py: "  << i_py
		    << " i_pz: "  << i_pz	      
		    << " re_linear 2D: " << i_px*p_bins[1] + i_py
		    << " re_linear 3D: " << (i_px*p_bins[1]*p_bins[2])+(i_py*p_bins[2])+i_pz	      
		    << " tot_w: " << tot_w
		    << " tot_px: " << tot_px 
		    << " tot_py: " << tot_py 
		    << " tot_pz: " << tot_pz 
		    << " rest_mo: " << mo
		    << std::endl; 
	}
	
	// Vranic et al. : epsilon_a, pa e2= (pc)2 + (m0c2)2
	double eps_a = tot_e/tot_w;
	double pa = std::sqrt(std::pow(eps_a,2)-std::pow(mo*C_LIGHT*C_LIGHT,2))/C_LIGHT;
	
	// Total p_norm
	double tot_p_norm = std::sqrt(
				      std::pow(tot_px,2)
				      +std::pow(tot_py,2)
				      +std::pow(tot_pz,2)
				      );
	
	// Vranic et al: angle between pa and pt, pb and pt 
	double cos_w = std::min(tot_p_norm / (tot_w*pa),1.0);
	double sin_w = std::sqrt(1 - cos_w*cos_w);
	
	// Inverse total p
	double inv_tot_p_norm = 1./tot_p_norm;
	
	// Computation of u1 unit vector
	double u1_x = tot_px*inv_tot_p_norm;
	double u1_y = tot_py*inv_tot_p_norm;
	double u1_z = tot_pz*inv_tot_p_norm; //0. in 2D case
	
	// Vranic et al. vec_d vector
	double d_vec_x = p_min[0] + (i_px+0.5)*dp[0];
	double d_vec_y = p_min[1] + (i_py+0.5)*dp[1];
	double d_vec_z = p_min[2] + (i_pz+0.5)*dp[2];//0. in 2D case
	
	// u3 = u1 x d_vec
	double u3_x = u1_y*d_vec_z - u1_z*d_vec_y; // 0. in 2D case 
	double u3_y = u1_z*d_vec_x - u1_x*d_vec_z; // 0. in 2D case
	double u3_z = u1_x*d_vec_y - u1_y*d_vec_x; // (!=0) ,along Z. dir
	
	
	  if ((m_verbosity>1) && (0 == m_mpi_rank)){	 
	    std::cout << " eps_a: " << eps_a
		      << " pa: " << pa
		      << " tot_p_norm: " << tot_p_norm 
		      << " cos_w: " << cos_w
		      << " sin_w: " << sin_w
	              << " u1_x: " << u1_x
	              << " u1_y: " << u1_y
	              << " u1_z: " << u1_z
	              << " d_x: " <<  d_vec_x
	              << " d_y: " <<  d_vec_y
	              << " d_z: " <<  d_vec_z
	              << " u3_x: " << u3_x
	              << " u3_y: " << u3_y
	              << " u3_z: " << u3_z	      
	              << std::endl; 
	  }
	  
	  // All particle momenta are not collinear
	  if (fabs(u3_x*u3_x + u3_y*u3_y + u3_z*u3_z) > 0)
	    {
	      
	      double u2_x = u1_y*u3_z - u1_z*u3_y;
	      double u2_y = u1_z*u3_x - u1_x*u3_z;
	      double u2_z = u1_x*u3_y - u1_y*u3_x;
	      
	      double u2_norm = 1./sqrt(u2_x*u2_x + u2_y*u2_y + u2_z*u2_z);
	      
	      // u2 normalized 
	      u2_x = u2_x * u2_norm;
	      u2_y = u2_y * u2_norm;
	      u2_z = u2_z * u2_norm;
	      
              // Select only first particle
	      // Tagging all others to be removed ...
	      
	      if ((m_verbosity>1) && (0 == m_mpi_rank)){
		std::cout << " cond. fabs() > 0 "
			  << " u2_x: " << u2_x
			  << " u2_y: " << u2_y
			  << " u2_z: " << u2_z
		          << std::endl; 
	      }
	      
	      // Update momentum of 2 first particles in the cell
              // First particle
	      pp.px[p_indexes[i_pcell][0]] = pa*(cos_w*u1_x + sin_w*u2_x);
	      pp.py[p_indexes[i_pcell][0]] = pa*(cos_w*u1_y + sin_w*u2_y);
	      pp.pz[p_indexes[i_pcell][0]] = pa*(cos_w*u1_z + sin_w*u2_z);
	      pp.w[p_indexes[i_pcell][0]]  = 0.5*tot_w;
	      
	      // Second particle
	      pp.px[p_indexes[i_pcell][1]] = pa*(cos_w*u1_x - sin_w*u2_x);
	      pp.py[p_indexes[i_pcell][1]] = pa*(cos_w*u1_y - sin_w*u2_y);
	      pp.pz[p_indexes[i_pcell][1]] = pa*(cos_w*u1_z - sin_w*u2_z);
	      pp.w[p_indexes[i_pcell][1]]  = 0.5*tot_w;
	      
	      // Mask the other indexes in the cell
	      for (int p=2; p<npart_per_cell; p++){
		m_mask_indexes.push_back(p_indexes[i_pcell][p]);
	      }
	    }//!(fabs()>0)
      }//!for(p-cells)
    }//! p_reduction
    
  }//! ns dim_3
}//! ns converter
