
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

#include "mpi.h"
#include "cxxopts.hpp"

#include <openPMD/openPMD.hpp>
using namespace openPMD;
using namespace std;

extern "C" {
  struct part{ int l_px, l_py, l_pz; double* px, *py, *pz;};    
  void read_particle(const char* cstring, int clen,
		     const char* spec,    int len_sp,
		     part* arrays, long* npart, long* npart_proc, long* start);
  void init_read();
}

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

void sdf_io(int argc, char *argv[]) {
  
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  std::string sdf_file;
  std::string species_name;
  
  cxxopts::Options optparse("sdf_io", "runs tests on sdf_io interface");
  optparse.add_options()(
			 "f,sdf_file", "sdf file to process",
                         cxxopts::value<std::string>(sdf_file))
                         ("s,species_name", "species_name",
                          cxxopts::value<std::string>(species_name));
  
  auto opts = optparse.parse(argc, argv);
  if (mpi_rank == 0 ){
    std::cout << "sdf2opmd_2d: file to process: " << sdf_file.c_str() << std::endl;
    std::cout << "sdf2opmd_2d: species name list: " << species_name.c_str() << std::endl;
  }
  
  // Get every species  
  std::vector<std::string> species_list = split(species_name.c_str());
  
  // Create the TTree
  std::string hdf_file=sdf_file.substr(0,sdf_file.find_last_of('.'));
  hdf_file +=".h5";
  
  // intialize Epoch read process 
  init_read();
  
  // Creating series
  Series series= Series(hdf_file.c_str(), Access::CREATE, MPI_COMM_WORLD);
  series.setAuthor("d.bertini@gsi.de");
  
  if (0 == mpi_rank)
    std::cout << " Creating Series with HDF output file: " << hdf_file << std::endl;
 
  // Loop for all species and fetch data in parallel
  for (size_t i=0; i < species_list.size(); i++) {
    std::string spec = species_list[i];
    if ( 0 == mpi_rank )
    std::cout<< " fetching data for  species: " << spec << std::endl; 

    part arrays;
    long e_npart{0}, e_npart_proc{0}, e_start{0};
    read_particle(sdf_file.c_str(), strlen(sdf_file.c_str()),
                  spec.c_str(),     strlen(spec.c_str()),  
		  &arrays, &e_npart, &e_npart_proc, &e_start);    
    
    if (0 == mpi_rank ) {
      // std::cout << "species: " << spec << " size: " << arrays.l_px << std::endl;
      //std::cout << "npart: " << e_npart << " npart_proc: " << e_npart_proc << " arrays_l: " << arrays.l_px
      //		<<   " start: " << e_start << std::endl;
      
      /*
      for (int k=0; k<arrays.l_px; k++){
	if ( ( (k%10000) == 0 ) && (mpi_rank == 0) )
	  std::cout << "rank: " << mpi_rank <<  " k: " << k << " px: " << arrays.px[k]
		    << " py: " << arrays.py[k] << " pz: " << arrays.pz[k] << std::endl; 
      }
      */
    }
    
    // Controlling Extension    
    assert( arrays.l_px == e_npart_proc );
    assert( arrays.l_py == e_npart_proc );
    assert( arrays.l_pz == e_npart_proc );
    
    // Create Particle species
    ParticleSpecies e = series.iterations[1].particles[spec.c_str()];
    
    // Create Dataset
    Datatype datatype = determineDatatype<double>();
    //Extent global_extent = {static_cast<unsigned long>(e_npart)};
    Extent global_extent = {e_npart};
    Dataset dataset = Dataset(datatype, global_extent);
    
    if (0 == mpi_rank)
      cout << "Prepared a Dataset of size " << dataset.extent[0] 
	   << " and Datatype " << dataset.dtype
	   << '\n';    
    
    Offset chunk_offset = {e_start-1};
    Extent chunk_extent = {e_npart_proc};
    
    e["momentum"]["x"].resetDataset(dataset);
    e["momentum"]["x"].storeChunk(shareRaw(arrays.px), chunk_offset, chunk_extent);

    e["momentum"]["y"].resetDataset(dataset);
    e["momentum"]["y"].storeChunk(shareRaw(arrays.py), chunk_offset, chunk_extent);
    
    e["momentum"]["z"].resetDataset(dataset);
    e["momentum"]["z"].storeChunk(shareRaw(arrays.pz), chunk_offset, chunk_extent);
    

    if (0 == mpi_rank)
      cout << "Registered a single chunk per MPI rank containing its "
	"contribution, "
	"ready to write content to disk\n";

    series.flush();

  } // species++

}

int main(int argc, char *argv[])
{
  int mpi_s{-1};
  int mpi_r{-1};
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_s);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_r);
  
  auto mpi_size = static_cast<uint64_t>(mpi_s);
  auto mpi_rank = static_cast<uint64_t>(mpi_r);
  
  sdf_io(argc, argv);
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}




    
