
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
  struct part{
    int   l_x, l_y, l_z, l_w, l_px, l_py, l_pz, l_vx, l_vy, l_vz, l_ek, l_rm, l_gm;
    double *x, *y, *z, *w,  *px, *py, *pz, *vx, *vy, *vz, *ek, *rm, *gm;
  };

  
  void read_particle(const char* cstring, int clen,
		     const char* spec,    int len_sp,
		     part* arrays, long* npart, long* npart_proc, long* start);

  struct field{
    int   global_sx, global_sy, global_sz, l_sx, l_sy, l_sz, l_dx, l_dy, l_dz, l_gridx, l_gridy, l_gridz, l_data_size;
    double stagger;
    double* gridx,  *gridy, *gridz, *l_field_data;
  };

  
  void read_fields( const char* cstring, int clen,  const char* blockid, int blen, field* field_x);
  void read_derived_vars( const char* cstring, int clen,  const char* blockid, int blen, field* field_x);  
  
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


int sdf_io(int argc, char *argv[]) {
  
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  std::string sdf_file;
  std::string info_name;
  std::string species_name;
  std::string field_name;
  std::string output_format;
  std::string derived_name;
  
  cxxopts::Options optparse("sdf_io", "runs tests on sdf_io interface");
  optparse.add_options()(
			 "f,sdf_file", "sdf file to process",
                         cxxopts::value<std::string>(sdf_file))
                         ("i,info_name", "info_name",
                         cxxopts::value<std::string>(info_name))    
                         ("m,field_name", "field_name",
                         cxxopts::value<std::string>(field_name))
                         ("d,derived_name", "derived_name",
                         cxxopts::value<std::string>(derived_name))    
                         ("s,species_name", "species_name",
                         cxxopts::value<std::string>(species_name))
                         ("o,output_format", "output_format",
                         cxxopts::value<std::string>(output_format));
  
  auto opts = optparse.parse(argc, argv);
  if (mpi_rank == 0 ){
    if ( info_name.length() == 0 ){
      std::cout << std::endl;
      std::cout << "sdf2opmd_3d: file to process: " << sdf_file.c_str() << std::endl;
      std::cout << "sdf2opmd_3d: field name list: " << field_name.c_str() << std::endl;
      std::cout << "sdf2opmd_3d: species name list: " << species_name.c_str() << std::endl;
      std::cout << "sdf2opmd_3d: derived name list: " << derived_name.c_str() << std::endl;    
      std::cout << "sdf2opmd_3d: output format: " << output_format.c_str() << std::endl;
      std::cout << std::endl;
    }else{
      std::cout << std::endl;
      std::cout << "sdf2opmd_3d: dry run:  reading SDF blocks info: " << info_name << " from file: " << sdf_file << std::endl;      
      std::cout << std::endl;
    }
  }

  // intialize Epoch read process 
  init_read();
  
  // First check if we want only SDF block info
  if (info_name.length() > 0 ){
    field field_x;
    read_derived_vars(sdf_file.c_str(), strlen(sdf_file.c_str()), info_name.c_str(), strlen(info_name.c_str()), &field_x);
    return 0;
  }

  // Get every field  
  std::vector<std::string> field_list;
  if ( !field_name.empty() )
    field_list= split(field_name.c_str());

  // Get every derived field  
  std::vector<std::string> derived_list;
  if ( !derived_name.empty() )
    derived_list = split(derived_name.c_str());
  
  // Get every species  
  std::vector<std::string> species_list;
  if ( !species_name.empty() )
    species_list = split(species_name.c_str());
  
  // Output filename
  std::string hdf_file=sdf_file.substr(0,sdf_file.find_last_of('.'));

  if ( output_format == "hdf5" ){
    hdf_file +=".h5";
  }else if ( output_format == "adios"){
    hdf_file +=".bp";
    // JSON not supported with MPI 
    //}else if ( output_format == "json"){
    // hdf_file +=".json";    
  } else {
    // default to hdf5
    hdf_file +=".h5";
  }

  if (0 == mpi_rank){
    std::cout << std::endl; 
    std::cout << "Creating Series with output file: " << hdf_file << std::endl;
  }
  
  // Creating series
  Series series= Series(hdf_file.c_str(), Access::CREATE, MPI_COMM_WORLD);
  series.setAuthor("d.bertini@gsi.de");
  
  // Read field data
  if (0 == mpi_rank){
    std::cout << "Reading field data from : " << sdf_file << std::endl;
    std::cout << std::endl;
  }
    
  //
  // Derived fields 
  //
  
  for (size_t k=0; k < derived_list.size(); k++) {
    std::string blockid = derived_list[k];
    
    if ( 0 == mpi_rank )
      std::cout<< " fetching derived data for field: " << blockid << std::endl; 
    
    field field_x;
    read_fields(sdf_file.c_str(), strlen(sdf_file.c_str()), blockid.c_str(), strlen(blockid.c_str()), &field_x);
    
    if (0 == mpi_rank){
      std::cout << std::endl;      
      std::cout << " C side: read_field: " << blockid.c_str() << std::endl;
      std::cout << " global sizes sx: " << field_x.global_sx << "  sy: " << field_x.global_sy << "  sz: " << field_x.global_sz << std::endl;
      std::cout << " local sizes l_sx: " << field_x.l_sx << " l_sy: " << field_x.l_sy << " l_sz: " << field_x.l_sz << std::endl;
      std::cout << " start sizes    st_x: " << field_x.l_dx << " st_y: " << field_x.l_dy << " st_z: " << field_x.l_dz << std::endl;
      std::cout << " grid  size  x: " << field_x.l_gridx << " y: " << field_x.l_gridy << " z: " << field_x.l_gridz << std::endl; 
      std::cout << std::endl;
    }

    // Loads the Grids  
    // Create a mesh for grids
    MeshRecordComponent gx_mesh =
      series.iterations[1].meshes["grid_x"][MeshRecordComponent::SCALAR];
    MeshRecordComponent gy_mesh =
      series.iterations[1].meshes["grid_y"][MeshRecordComponent::SCALAR];
    MeshRecordComponent gz_mesh =
      series.iterations[1].meshes["grid_z"][MeshRecordComponent::SCALAR];
        
    Datatype dtype = determineDatatype<double>();;
    Extent extent_x = {field_x.l_gridx};
    Dataset dataset_x = Dataset(dtype, extent_x);
    Extent extent_y = {field_x.l_gridy};
    Dataset dataset_y = Dataset(dtype, extent_y);
    Extent extent_z = {field_x.l_gridz};
    Dataset dataset_z = Dataset(dtype, extent_z);
    
    
    gx_mesh.resetDataset(dataset_x);
    gy_mesh.resetDataset(dataset_y);
    gz_mesh.resetDataset(dataset_z);    
    
    Offset offset_x = {0};
    Offset offset_y = {0};
    Offset offset_z = {0};
    
    gx_mesh.storeChunkRaw( (field_x.gridx), offset_x, extent_x);
    gy_mesh.storeChunkRaw( (field_x.gridy), offset_y, extent_y);      
    gz_mesh.storeChunkRaw( (field_x.gridz), offset_z, extent_z);
    
    // Load fields values 
    // Remap in 2D using vector
    const int nx = field_x.l_sx;
    const int ny = field_x.l_sy;
    const int nz = field_x.l_sz;
    std::vector<double> v_map;
    
    // Fill the 3D vector map
    for (int i=0;i<nx; i++){
      for (int j=0;j<ny; j++){
	for (int kk=0;kk<nz; kk++){
	  int l_index = (kk*nx*ny) + (j*nx) + (i);
	  v_map.push_back(field_x.l_field_data[l_index]);
	}
      }
    }
    
    // print the 2D map
    //for(int i = 0; i < nx*ny; i++) std::cout << " i: " << i << " val: " <<  v_map[i] << std::endl; 

    std::replace(blockid.begin(), blockid.end(), '/', '_');
    
    // Create a 2D mesh for E_x
    MeshRecordComponent ex_mesh =
      series.iterations[1].meshes[blockid.c_str()][MeshRecordComponent::SCALAR];
    
    // Only 1D domain decomposition in first index
    Datatype datatype = determineDatatype<double>();
    Extent global_extent = {field_x.global_sx, field_x.global_sy, field_x.global_sz};
    Dataset dataset = Dataset(datatype, global_extent);
    
    
    if (0 == mpi_rank)
      cout << "Prepared a Field Dataset of size " << dataset.extent[0] << "x"
	   << dataset.extent[1] << " and Datatype " << dataset.dtype
	   << '\n';
    
    ex_mesh.resetDataset(dataset);
    
    // example shows a 1D domain decomposition in first index
    Offset chunk_offset = {field_x.l_dx, field_x.l_dy, field_x.l_dz};
    Extent chunk_extent = {field_x.l_sx, field_x.l_sy, field_x.l_sz};
    
    
    // Store the 2D map
    ex_mesh.storeChunk(v_map, chunk_offset, chunk_extent);
    
    if (0 == mpi_rank)
      cout << "Registered a single chunk per MPI rank containing its "
	"contribution, "
	"ready to write content to disk\n";
    
    series.flush();
    
    if (0 == mpi_rank)
      cout << "Fields Datasets content has been fully written to disk\n";
  } // fields++

  
  //
  // Fields
  //
    
  for (size_t k=0; k < field_list.size(); k++) {
    std::string blockid = field_list[k];
    if ( 0 == mpi_rank )
      std::cout<< " fetching data for field: " << blockid << std::endl; 
    
    field field_x;
    read_fields(sdf_file.c_str(), strlen(sdf_file.c_str()), blockid.c_str(), strlen(blockid.c_str()), &field_x);
    
    if (0 == mpi_rank){
      std::cout << std::endl;      
      std::cout << " C side: read_field: " << blockid.c_str() << std::endl;
      std::cout << " global sizes sx: " << field_x.global_sx << "  sy: " << field_x.global_sy << "  sz: " << field_x.global_sz << std::endl;
      std::cout << " local sizes l_sx: " << field_x.l_sx << " l_sy: " << field_x.l_sy << " l_sz: " << field_x.l_sz << std::endl;
      std::cout << " start sizes    st_x: " << field_x.l_dx << " st_y: " << field_x.l_dy << " st_z: " << field_x.l_dz << std::endl;
      std::cout << " grid  size  x: " << field_x.l_gridx << " y: " << field_x.l_gridy << " z: " << field_x.l_gridz << std::endl; 
      std::cout << std::endl;
    }

    // Loads the Grids  
    // Create a mesh for grids
    MeshRecordComponent gx_mesh =
      series.iterations[1].meshes["grid_x"][MeshRecordComponent::SCALAR];
    MeshRecordComponent gy_mesh =
      series.iterations[1].meshes["grid_y"][MeshRecordComponent::SCALAR];
    MeshRecordComponent gz_mesh =
      series.iterations[1].meshes["grid_z"][MeshRecordComponent::SCALAR];    
        
    Datatype dtype = determineDatatype<double>();;
    Extent extent_x = {field_x.l_gridx};
    Dataset dataset_x = Dataset(dtype, extent_x);
    Extent extent_y = {field_x.l_gridy};
    Dataset dataset_y = Dataset(dtype, extent_y);
    Extent extent_z = {field_x.l_gridz};
    Dataset dataset_z = Dataset(dtype, extent_z);    
    
    gx_mesh.resetDataset(dataset_x);
    gy_mesh.resetDataset(dataset_y);
    gz_mesh.resetDataset(dataset_z);
    
    Offset offset_x = {0};
    Offset offset_y = {0};
    Offset offset_z = {0};    
    
    gx_mesh.storeChunkRaw( (field_x.gridx), offset_x, extent_x);
    gy_mesh.storeChunkRaw( (field_x.gridy), offset_y, extent_y);
    gz_mesh.storeChunkRaw( (field_x.gridz), offset_z, extent_z);          
    
    // Load fields values 
    // Remap in 2D using vector
    const int nx = field_x.l_sx;
    const int ny = field_x.l_sy;
    const int nz = field_x.l_sz;    
    std::vector<double> v_map;
    
    // Fill the 2D vector map
    for (int i=0;i<nx; i++){
      for (int j=0;j<ny; j++){
	for (int kk=0;kk<nz; kk++){
	  int l_index = (kk*nx*ny) + (j*nx) + (i);
	  v_map.push_back(field_x.l_field_data[l_index]);
	}
      }
    }

    
    // print the 2D map
    //for(int i = 0; i < nx*ny; i++) std::cout << " i: " << i << " val: " <<  v_map[i] << std::endl; 
    
    // Create a 2D mesh for E_x
    MeshRecordComponent ex_mesh =
      series.iterations[1].meshes[blockid.c_str()][MeshRecordComponent::SCALAR];
    
    // Only 1D domain decomposition in first index
    Datatype datatype = determineDatatype<double>();
    Extent global_extent = {field_x.global_sx, field_x.global_sy, field_x.global_sz};
    Dataset dataset = Dataset(datatype, global_extent);
    
    
    if (0 == mpi_rank)
      cout << "Prepared a Field Dataset of size " << dataset.extent[0] << "x"
	   << dataset.extent[1] << " and Datatype " << dataset.dtype
	   << '\n';
    
    ex_mesh.resetDataset(dataset);
    
    // example shows a 1D domain decomposition in first index
    Offset chunk_offset = {field_x.l_dx, field_x.l_dy, field_x.l_dz};
    Extent chunk_extent = {field_x.l_sx, field_x.l_sy, field_x.l_sz};
    
    
    // Store the 2D map
    ex_mesh.storeChunk(v_map, chunk_offset, chunk_extent);
    
    if (0 == mpi_rank)
      cout << "Registered a single chunk per MPI rank containing its "
	"contribution, "
	"ready to write content to disk\n";
    
    series.flush();
    
    if (0 == mpi_rank)
      cout << "Fields Datasets content has been fully written to disk\n";
  } // fields++


  // 
  // Particles
  //
  
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

    /*
    if (0 == mpi_rank ) {
      std::cout << "species: " << spec << " size: " << arrays.l_px << std::endl;
      std::cout << "npart: " << e_npart << " npart_proc: " << e_npart_proc << " arrays_l: " << arrays.l_px
      		<<   " start: " << e_start << std::endl;
      
      for (int k=0; k<arrays.l_px; k++){
	if ( ( (k%10000) == 0 ) && (mpi_rank == 0) )
	  std::cout << "rank: " << mpi_rank <<  " k: " << k << " px: " << arrays.px[k]
		    << " py: " << arrays.py[k] << " pz: " << arrays.pz[k] << std::endl; 
      }    
    }
    */    

    // Controlling Extensions
    if (arrays.l_px>0) 
    assert( arrays.l_px == e_npart_proc );
    if (arrays.l_py>0) 
    assert( arrays.l_py == e_npart_proc );
    if (arrays.l_pz>0) 
    assert( arrays.l_pz == e_npart_proc );
    if (arrays.l_x>0) 
    assert( arrays.l_x  == e_npart_proc );
    if (arrays.l_y>0) 
    assert( arrays.l_y  == e_npart_proc );
    if (arrays.l_z>0) 
    assert( arrays.l_z  == e_npart_proc );    
    if (arrays.l_w>0) 
    assert( arrays.l_w  == e_npart_proc );
    
    
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

    // Add particules infos accroding to input SDF
    
    // P
    if (arrays.l_px>0){ 
      e["momentum"]["x"].resetDataset(dataset);
      e["momentum"]["x"].storeChunkRaw( (arrays.px), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_py>0){ 
      e["momentum"]["y"].resetDataset(dataset);
      e["momentum"]["y"].storeChunkRaw( (arrays.py), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_pz>0){     
      e["momentum"]["z"].resetDataset(dataset);
      e["momentum"]["z"].storeChunkRaw( (arrays.pz), chunk_offset, chunk_extent);
    }
    
    // V
    if (arrays.l_vx>0){     
      e["velocity"]["x"].resetDataset(dataset);
      e["velocity"]["x"].storeChunkRaw( (arrays.vx), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_vy>0){     
      e["velocity"]["y"].resetDataset(dataset);
      e["velocity"]["y"].storeChunkRaw( (arrays.vy), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_vz>0){         
      e["velocity"]["z"].resetDataset(dataset);
      e["velocity"]["z"].storeChunkRaw( (arrays.vz), chunk_offset, chunk_extent);
    }
    
    // X
    if (arrays.l_x>0){         
      e["position"]["x"].resetDataset(dataset);
      e["position"]["x"].storeChunkRaw( (arrays.x), chunk_offset, chunk_extent);
    }
    
    //Y
    if (arrays.l_y>0){         
      e["position"]["y"].resetDataset(dataset);
      e["position"]["y"].storeChunkRaw( (arrays.y), chunk_offset, chunk_extent);
    }

    //Z
    if (arrays.l_z>0){         
      e["position"]["z"].resetDataset(dataset);
      e["position"]["z"].storeChunkRaw( (arrays.z), chunk_offset, chunk_extent);
    }

    auto const scalar = openPMD::RecordComponent::SCALAR;
    
    if (arrays.l_w>0){         
      e["weighting"][scalar].resetDataset(dataset);
      e["weighting"][scalar].storeChunkRaw ( (arrays.w), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_ek>0){             
      e["ek"][scalar].resetDataset(dataset);
      e["ek"][scalar].storeChunkRaw( (arrays.ek), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_rm>0){             
      e["rm"][scalar].resetDataset(dataset);
      e["rm"][scalar].storeChunkRaw( (arrays.rm), chunk_offset, chunk_extent);
    }
    
    if (arrays.l_gm>0){             
      e["gm"][scalar].resetDataset(dataset);
      e["gm"][scalar].storeChunkRaw( (arrays.gm), chunk_offset, chunk_extent);    
    }  
    
    if (0 == mpi_rank)
      cout << "Registered a single chunk per MPI rank containing its "
	"contribution, "
	"ready to write content to disk\n";

    series.flush();

  } // species++

  return 0;
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

  if ( 0 == mpi_rank ) {
    std::cout << std::endl;
    std::cout <<"sdf2opmd_3d:: MPI initialized with size: " << mpi_size << " : " << mpi_rank << std::endl;
    std::cout << std::endl;
  }

  sdf_io(argc, argv);
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}




    
