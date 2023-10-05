#pragma once

#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <sys/types.h>
#include <thread>
#include <unistd.h>
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>

// MPI
#include "mpi.h"

// Options
#include "cxxopts.hpp"

// openPMD
#include <openPMD/openPMD.hpp>

using namespace openPMD;
using namespace std;

namespace converter
{
  namespace dim_2
  {  
    
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
    
    namespace si_units{
      enum si_units_t
	{
	  length = 0, // L
	  mass = 1, // M
	  time = 2, // T
	  electricCurrent = 3, // I
	  thermodynamicTemperature = 4, // theta
	  amountOfSubstance = 5, // N
	  luminousIntensity = 6, // J
	};
    }
    
    inline
    std::vector<std::string> split(const char *str, char c = ':')
    {
      std::vector<std::string> result;
      do
	{
	  const char *begin = str;
	  while(*str != c && *str)
	    str++;
	  result.push_back(std::string(begin, str));
	} while (0 != *str++);
      
      return result;
    }
    
    // Fields blockid value map
    enum class field_value { ex,
			     ey,
			     ez,
			     bx,
			     by,
			     bz,
			     jx,
			     jy,
			     jz			
    };
    
    static std::map<std::string, field_value > m_field_value
      {
	{"ex", field_value::ex }
	,{"ey", field_value::ey }
	,{"ez", field_value::ez }
	,{"bx", field_value::bx }
	,{"by", field_value::by }
	,{"bz", field_value::bz }
	,{"jx", field_value::jx }
	,{"jy", field_value::jy }
	,{"jz", field_value::jz }        
      };
    
    
    class Converter2D{  
      
    public:
      int sdf_io(int argc, char *argv[]);  
      
      inline
      void storeFieldRecord( openPMD::RecordComponent& field_mesh,
			     std::vector<double>& vmap,
			     const field& vfield
			     )
      {
	// Store as double     
	Datatype datatype = determineDatatype<double>();
	Extent global_extent = {vfield.global_sx, vfield.global_sy};
	Dataset dataset = Dataset(datatype, global_extent);
	
	field_mesh.resetDataset(dataset);
	
	// Chunk definition : Offset + Extent
	Offset chunk_offset = {vfield.l_dx, vfield.l_dy};
	Extent chunk_extent = {vfield.l_sx, vfield.l_sy};
	
	// Store the 2D map
	field_mesh.storeChunk(vmap, chunk_offset, chunk_extent);
      }
      
      bool isE(field_value val){
	if ( 
	    ( val == field_value::ex )
	    ||
	    ( val == field_value::ey )
	    ||
	    ( val == field_value::ez )
	     ) return true;
	
	return false;  
      }
      
      bool isB(field_value val){
	if ( 
	    ( val == field_value::bx )
	    ||
	    ( val == field_value::by )
	    ||
	    ( val == field_value::bz )
	     ) return true;
	
	return false;  
      }
      
      bool isJ(field_value val){
	if ( 
	    ( val == field_value::jx )
	    ||
	    ( val == field_value::jy )
	    ||
	    ( val == field_value::jz )
	     ) return true;
	
	return false;  
      }
      
      struct fieldE{
      public:
	static std::vector<float> getUnitDimension()
	{
	  /* E is in volts per meters: V / m = kg * m / (A * s^3)
	   *   -> L * M * T^-3 * I^-1
	   */
	  std::vector<float> unitDimension(7, 0.0);
	  unitDimension.at(si_units::length) = 1.0;
	  unitDimension.at(si_units::mass) = 1.0;
	  unitDimension.at(si_units::time) = -3.0;
	  unitDimension.at(si_units::electricCurrent) = -1.0;
	  return unitDimension;
	}
      };
      
      struct fieldB{
      public:
	static std::vector<float> getUnitDimension()
	{
	  /* L, M, T, I, theta, N, J
	   *
	   * B is in Tesla : kg / (A * s^2)
	   *   -> M * T^-2 * I^-1
	   */
	  std::vector<float> unitDimension( 7, 0.0 );
	  unitDimension.at(si_units::mass) =  1.0;
	  unitDimension.at(si_units::time) = -2.0;
	  unitDimension.at(si_units::electricCurrent) = -1.0;	
	  return unitDimension;
	}
      };
      
      struct fieldJ{
      public:
	static std::vector<float> getUnitDimension()
	{
	  /* L, M, T, I, theta, N, J
	   *
	   * J is in A/m^2
	   *   -> L^-2 * I
	   */
	  std::vector<float> unitDimension( 7, 0.0 );
	  unitDimension.at(si_units::length) = -2.0;
	  unitDimension.at(si_units::electricCurrent) =  1.0;
	  return unitDimension;	
	}
      };     
      
      
      inline void
      setOpmdUnits ( openPMD::Mesh mesh, const std::string field_name )
      {
	if (field_name[0] == 'E'){  // Electric field
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::L,  1},
	      {openPMD::UnitDimension::M,  1},
	      {openPMD::UnitDimension::T, -3},
	      {openPMD::UnitDimension::I, -1},
	    });
	} else if (field_name[0] == 'B'){ // Magnetic field
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::M,  1},
	      {openPMD::UnitDimension::I, -1},
	      {openPMD::UnitDimension::T, -2}
	    });
	} else if (field_name[0] == 'j'){ // current
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::L, -2},
	      {openPMD::UnitDimension::I,  1},
	    });
	} else if (field_name.substr(0,3) == "rho"){ // charge density
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::L, -3},
	      {openPMD::UnitDimension::I,  1},
	      {openPMD::UnitDimension::T,  1},
	    });
	}
      }
      
      inline std::map< openPMD::UnitDimension, double >
      getOpmdUnits ( std::string const & record_name )
      {
	if( (record_name == "position") || (record_name == "positionOffset") ) return {
	    {openPMD::UnitDimension::L,  1.}
	  };
	else if( record_name == "momentum" ) return {
	    {openPMD::UnitDimension::L,  1.},
	    {openPMD::UnitDimension::M,  1.},
	    {openPMD::UnitDimension::T, -1.}
	  };
	else if( record_name == "charge" ) return {
	    {openPMD::UnitDimension::T,  1.},
	    {openPMD::UnitDimension::I,  1.}
	  };
	else if( record_name == "mass" ) return {
	    {openPMD::UnitDimension::M,  1.}
	  };
	else if( record_name == "E" ) return {
	    {openPMD::UnitDimension::L,  1.},
	    {openPMD::UnitDimension::M,  1.},
	    {openPMD::UnitDimension::T, -3.},
	    {openPMD::UnitDimension::I, -1.},
	  };
	else if( record_name == "B" ) return {
	    {openPMD::UnitDimension::M,  1.},
	    {openPMD::UnitDimension::I, -1.},
	    {openPMD::UnitDimension::T, -2.}
	  };	  
	else if( record_name == "J" ) return {
	    {openPMD::UnitDimension::L, -2},
	    {openPMD::UnitDimension::I,  1}
	  };
	else if( record_name == "rho" ) return {
	    {openPMD::UnitDimension::L, -3},
	    {openPMD::UnitDimension::I,  1},
	    {openPMD::UnitDimension::T,  1},
	  };	    
	
	else return {}; 
      }
      
    }; // ! class Converter2D
    
  } //ns_2d
}// ns_converter  
