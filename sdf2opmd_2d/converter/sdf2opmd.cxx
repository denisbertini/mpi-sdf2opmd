
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

extern "C" {
  struct part{ int l_px, l_py, l_pz; double* px, *py, *pz;};    
  void read_particle(const char* cstring, int clen,
		     const char* spec,    int len_sp,
		     part* arrays);
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
    
  int rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
  if (rank == 0 ){
    std::cout << "sdf2opmd_2d: file to process: " << sdf_file.c_str() << std::endl;
    std::cout << "sdf2opmd_2d: species name list: " << species_name.c_str() << std::endl;
  }

  // Get every species  
  std::vector<std::string> species_list = split(species_name.c_str());

  // Create the TTree
  std::string hdf_file=sdf_file.substr(0,sdf_file.find_last_of('.'));
  //root_file += "_r"+std::to_string(rank);
  hdf_file +=".h5";

  std::cout << " Creating Series with HDF output file: " << hdf_file << std::endl;
  
  // open file for writing
 
  Series series = Series(
  			 hdf_file.c_str(), Access::CREATE, MPI_COMM_WORLD);

  std::cout << "Series  created with HDF output file: " << hdf_file << std::endl;
 
  // Initialise ONCE the setup for MPI-I/O
  init_read();

  // Loop for all species and fetch data in parallel
  for (size_t i=0; i < species_list.size(); i++) {
    std::string spec = species_list[i];
    std::cout<< " fetching data for  species: " << spec << std::endl; 

    // Particle Kinematics
    double px[species_list.size()];
    double py[species_list.size()];
    double pz[species_list.size()];
    
    part arrays;
    read_particle(sdf_file.c_str(), strlen(sdf_file.c_str()),
                  spec.c_str(),     strlen(spec.c_str()),  
		  &arrays);    
    
    //std::cout << "species: " << spec << " size: " << arrays.l_px << std::endl;
    
    for (int k=0; k<arrays.l_px; k++){
      if ( ( (k%10000) == 0 ) && (rank == 0) )
	std::cout << "rank: " << rank <<  " k: " << k << " px: " << arrays.px[k]
		  << " py: " << arrays.py[k] << " pz: " << arrays.pz[k] << std::endl; 
      px[i] = arrays.px[k];
      py[i] = arrays.py[k];
      pz[i] = arrays.pz[k];      
    }

    /*
    
    // Now copy back to vector ? do i need that?
    std::vector<double> v_px(px, px + species_list.size());
    std::vector<double> v_py(py, py + species_list.size());
    std::vector<double> v_pz(pz, pz + species_list.size());

    ParticleSpecies part = series.iterations[0].particles[spec.c_str()];

    Datatype dtype = determineDatatype(shareRaw(v_px));
    Extent size = {mpi_size * v_px.size()};
    auto dataset = Dataset(dtype, size, "{ \"resizable\": false }");

    RecordComponent rc_x = part["momentum"]["x"];
    RecordComponent rc_y = part["momentum"]["y"];
    RecordComponent rc_z = part["momentum"]["z"];
 
    rc_x.resetDataset(dataset);
    rc_y.resetDataset(dataset);
    rc_z.resetDataset(dataset);    

    Offset offset = {0};
    rc_x.storeChunk(v_px, offset, {v_px.size()});
    rc_y.storeChunk(v_py, offset, {v_py.size()});
    rc_x.storeChunk(v_pz, offset, {v_pz.size()});

    // do I/O
    series.flush();
    if ( rank == 0 )
    std::cout << "Dataset written to output HDF file: "<< hdf_file << std::endl;
    
    */
  }//!species

 std::cout << "exiting sdfio(): "  << std::endl;

}

int main(int argc, char *argv[]) {

  auto start = std::chrono::high_resolution_clock::now();

  int mpi_rank, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    {
        // global data set to write: [MPI_Size * 10, 300]
        // each rank writes a 10x300 slice with its MPI rank as values
        auto const value = float(mpi_size);
        std::vector<float> local_data(10 * 300, value);
        if (0 == mpi_rank)
	  std::cout << "Set up a 2D array with 10x300 elements per MPI rank ("
                 << mpi_size << "x) that will be written to disk\n";

        // open file for writing
        Series series = Series(
            "/lustre/rz/dbertini/parallel_write.h5", Access::CREATE, MPI_COMM_WORLD);
        if (0 == mpi_rank)
	  std::cout << "Created an empty series in parallel with " << mpi_size
                 << " MPI ranks\n";

        MeshRecordComponent mymesh =
            series.iterations[1].meshes["mymesh"][MeshRecordComponent::SCALAR];

        // example 1D domain decomposition in first index
        Datatype datatype = determineDatatype<float>();
        Extent global_extent = {10ul * mpi_size, 300};
        Dataset dataset = Dataset(datatype, global_extent);

        if (0 == mpi_rank)
	  std::cout << "Prepared a Dataset of size " << dataset.extent[0] << "x"
                 << dataset.extent[1] << " and Datatype " << dataset.dtype
                 << '\n';

        mymesh.resetDataset(dataset);
        if (0 == mpi_rank)
	  std::cout << "Set the global Dataset properties for the scalar field "
                    "mymesh in iteration 1\n";

        // example shows a 1D domain decomposition in first index
        Offset chunk_offset = {10ul * mpi_rank, 0};
        Extent chunk_extent = {10, 300};
        mymesh.storeChunk(local_data, chunk_offset, chunk_extent);
        if (0 == mpi_rank)
	  std::cout << "Registered a single chunk per MPI rank containing its "
                    "contribution, "
                    "ready to write content to disk\n";

        series.flush();
        if (0 == mpi_rank)
	  std::cout << "Dataset content has been fully written to disk\n";
    }

    // openPMD::Series MUST be destructed at this point
    MPI_Finalize();

    return 0;

}

