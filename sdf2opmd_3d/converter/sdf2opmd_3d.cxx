#include <map>
#include <algorithm>
#include "sdf2opmd_3d.h"
#include "merger_3d.h"

namespace converter
{
  namespace dim_3
  {  

    namespace fs = std::filesystem;
    
    int Converter3D::sdf_io(int argc, char *argv[]) {
      int mpi_rank, mpi_size;
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      
      std::string sdf_dir;
      std::string sdf_files;
      std::string info_name;
      std::string species_name;
      std::string compression_name;      
      std::string field_name;
      std::string output_format;
      std::string derived_name;
      
      cxxopts::Options optparse("sdf_io", "runs tests on sdf_io interface");
      optparse.add_options()(
	 "p,dir", "sdf directory to process",
         cxxopts::value<std::string>(sdf_dir))
	("f,sdf_file", "sdf file to process",
	 cxxopts::value<std::string>(sdf_files))
	("i,info_name", "info_name",
	 cxxopts::value<std::string>(info_name))    
	("m,field_name", "field_name",
	 cxxopts::value<std::string>(field_name))
	("d,derived_name", "derived_name",
	 cxxopts::value<std::string>(derived_name))    
	("s,species_name", "species_name",
	 cxxopts::value<std::string>(species_name))
	("z,compression", "compression",
	 cxxopts::value<std::string>(compression_name))	
	("o,output_format", "output_format",
	 cxxopts::value<std::string>(output_format));
      
      auto opts = optparse.parse(argc, argv);
      
      if (mpi_rank == 0 ){
	if ( info_name.length() == 0 ){
	  std::cout << std::endl;
	  std::cout << "sdf2opmd_3d: directory to process: " << sdf_dir.c_str() << std::endl;
	  std::cout << "sdf2opmd_3d: files to process: " << sdf_files.c_str() << std::endl;
	  std::cout << "sdf2opmd_3d: field name list: " << field_name.c_str() << std::endl;
	  std::cout << "sdf2opmd_3d: species name list: " << species_name.c_str() << std::endl;
	  std::cout << "sdf2opmd_3d: derived name list: " << derived_name.c_str() << std::endl;
	  std::cout << "sdf2opmd_3d: compression: " << compression_name.c_str() << std::endl;	  
	  std::cout << "sdf2opmd_3d: output format: " << output_format.c_str() << std::endl;
	  std::cout << std::endl;
	}else{
	  std::cout << std::endl;
	  std::cout << "sdf2opmd_3d: dry run:  reading SDF blocks info: " << info_name
		    << " from files: " << sdf_files << std::endl;      
	  std::cout << std::endl;
	}
      }
      
      // Initialize Epoch read process 
      init_read();

      // Get every files 
      std::vector<std::string> sdf_file_list;
      if ( !sdf_files.empty() ){
	sdf_file_list= split(sdf_files.c_str());
        for (unsigned int i=0;i<sdf_file_list.size();i++){
	  sdf_file_list[i]=sdf_dir+"/"+ sdf_file_list[i];
	}
      }
	
      std::vector<std::string> common_sdf_files = get_common_files(sdf_file_list, sdf_dir);

      // Main loop over sdf files 
      for(std::string sdf_file : common_sdf_files)
	{
	  if (0 == mpi_rank)
	    std::cout << "Converting SDF file: " << sdf_file << " to openPMD format ..." << std::endl;      
	  
	  // First check if we want only SDF block info
	  if (info_name.length() > 0 ){
	    field field_x;
	    read_derived_vars(sdf_file.c_str(), strlen(sdf_file.c_str()),
			      info_name.c_str(), strlen(info_name.c_str()), &field_x);
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
	  std::string digits=hdf_file.substr(hdf_file.size() - 4);
	  int sdf_iter=stoi(digits);      
	  
	  // Converted files dumped to local
	  if ( output_format == "hdf5" ){
	    hdf_file="%04T.h5";
	  }else if ( output_format == "adios"){
	    hdf_file="%04T.bp";
	  } else {
	    // default to hdf5
	    hdf_file="%04T.h5";	
	  }

	  if (0 == mpi_rank){
	    std::cout << std::endl; 
	    std::cout << "Creating Series with output file: " << hdf_file << std::endl;
	  }
	  
          //
	  // Open_pmd I/O
	  //
	  
	  // Creating series
	  Series series= Series(hdf_file.c_str(), Access::CREATE, MPI_COMM_WORLD);
	  series.setAuthor("d.bertini@gsi.de");
	  series.setMachine("Virgo2");
	  series.setSoftwareDependencies("https://git.gsi.de/d.bertini/pp-containers:prod/rlx8_ompi_ucx.def");
	  series.setMeshesPath("fields/");
	  series.setParticlesPath("particles/");
	  series.setIterationEncoding(IterationEncoding::fileBased);
	  series.setIterationFormat("%04T.h5");
	  
	  // In parallel contexts, it's important to explicitly open iterations.
	  series.iterations[sdf_iter].open();
	  
	  // Read field data
	  if (0 == mpi_rank){
	    std::cout << "Reading field data from : " << sdf_file << std::endl;
	    std::cout << std::endl;
	  }
	  
	  //
	  // Fields
	  //
	  
	  bool load_grids  = true;
	  bool meta_fields[3] = {true, true, true};
	  
	  for (size_t k=0; k < field_list.size(); k++) {
	    std::string blockid = field_list[k];
	    if ( 0 == mpi_rank )
	      std::cout<< "Fetching data for field: " << blockid << std::endl; 
	    
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
	    
	    // Create Mesh for Grids once 
	    if ( load_grids ) {  
	      if ( 0 == mpi_rank )
		std::cout<< " converting domain grids ..."<< std::endl; 
	      
	      MeshRecordComponent gx_mesh =
		series.iterations[sdf_iter].meshes["grid_x"][MeshRecordComponent::SCALAR];
	      MeshRecordComponent gy_mesh =
		series.iterations[sdf_iter].meshes["grid_y"][MeshRecordComponent::SCALAR];
	      MeshRecordComponent gz_mesh =
		series.iterations[sdf_iter].meshes["grid_z"][MeshRecordComponent::SCALAR];	  
	      
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
	      
	      load_grids=false;
	    }
	    
	    //
	    // Load fields values 
	    //
	    
	    // Remap in 3D using vector
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
	    
	    // print the 3D map
	    //for(int i = 0; i < nx*ny; i++) std::cout << " i: " << i << " val: " <<  v_map[i] << std::endl; 
	    
	    // Create a 3D mesh for Fields E, B, J 
	    auto meshes = series.iterations[sdf_iter].meshes;   
	    
	    bool flag  = false;
	    flag  = isE(m_field_value[blockid]);
	    
	    if (flag && meta_fields[0]){
	      auto mesh = meshes["E"];
	      mesh.setAxisLabels({"x", "y", "z"});
	      setOpmdUnits(mesh, "E");
	      double gx_spacing = (field_x.gridx[field_x.l_gridx-1]-field_x.gridx[0])/field_x.global_sx;
	      double gy_spacing = (field_x.gridy[field_x.l_gridy-1]-field_x.gridy[0])/field_x.global_sy;
	      double gz_spacing = (field_x.gridz[field_x.l_gridz-1]-field_x.gridz[0])/field_x.global_sz;
	      
	      mesh.setGridSpacing(std::vector<double>{gx_spacing, gy_spacing, gz_spacing})
		.setGridGlobalOffset({field_x.gridx[0], field_x.gridy[0], field_x.gridz[0]});	  
	      // Fields meta data loaded once
	      meta_fields[0]=false;
	    }
	    
	    flag  = isB(m_field_value[blockid]);
	    if (flag && meta_fields[1]){
	      auto mesh = meshes["B"];
	      mesh.setAxisLabels({"x", "y", "z"});
	      meta_fields[1]=false;
	    }
	    
	    flag  = isJ(m_field_value[blockid]);
	    if (flag && meta_fields[2]){
	      auto mesh = meshes["J"];
	      mesh.setAxisLabels({"x", "y", "z"});
	      meta_fields[1]=false;
	    }     
	    
	    if ( m_field_value[blockid] == field_value::ex ){       
	      auto mesh_ex = meshes["E"]["x"];
	      storeFieldRecord( mesh_ex, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::ey ) {
	      auto mesh_ey = meshes["E"]["y"];
	      storeFieldRecord( mesh_ey, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::ez ) {
	      auto mesh_ez = meshes["E"]["z"];
	      storeFieldRecord( mesh_ez, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::bx ) {
	      auto mesh_bx = meshes["B"]["x"];
	      storeFieldRecord( mesh_bx, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::by ) {
	      auto mesh_by = meshes["B"]["y"];
	      storeFieldRecord( mesh_by, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::bz ) {
	      auto mesh_bz = meshes["B"]["z"];
	      storeFieldRecord( mesh_bz, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::jx ) {
	      auto mesh_jx = meshes["J"]["x"];
	      storeFieldRecord( mesh_jx, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::jy ) {
	      auto mesh_jy = meshes["J"]["y"];
	      storeFieldRecord( mesh_jy, v_map, field_x);          
	    } else if ( m_field_value[blockid] == field_value::jz ) {
	      auto mesh_jz = meshes["J"]["z"];
	      storeFieldRecord( mesh_jz, v_map, field_x);          
	    }              
	    
	    // Sync to Disk    
	    series.flush();
	    
	    if (0 == mpi_rank)
	      cout << "Fields Datasets content has been fully written to disk\n";
	  } // fields++
	  
	  
	  //
	  // Derived fields 
	  //
	  
	  for (size_t k=0; k < derived_list.size(); k++) {
	    std::string blockid = derived_list[k];
	    
	    if ( 0 == mpi_rank )
	      std::cout<< "Fetching derived data for field: " << blockid << std::endl; 
	    
	    field field_x;
	    read_fields(sdf_file.c_str(), strlen(sdf_file.c_str()), blockid.c_str(), strlen(blockid.c_str()), &field_x);
	    
	    // Load fields values 
	    // Remap in 3D using vector
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
	    // Modify field semantics for file layout
	    std::replace(blockid.begin(), blockid.end(), '/', '_');
	    
	    // Create a 2D mesh for derived Fields  
	    auto meshes = series.iterations[sdf_iter].meshes;
	    auto mesh_d = meshes["derived"][blockid.c_str()];
	    storeFieldRecord( mesh_d, v_map, field_x);   
	    
	    // Sync to Disk
	    series.flush();
	    
	    if (0 == mpi_rank)
	      cout << "Derived Fields Datasets content has been fully written to disk\n";
	  } // derived fields++


	  // 
	  // Particles
	  //
	  
	  // Loop for all species and fetch data in parallel
	  for (size_t i_spec=0; i_spec < species_list.size(); i_spec++) {
	    std::string spec = species_list[i_spec];

	    if ( 0 == mpi_rank )
	      std::cout<< " fetching data for  species: " << spec << std::endl; 
	    
	    part arrays;
	    long e_npart{0}, e_npart_proc{0}, e_start{0};
	    read_particle(sdf_file.c_str(), strlen(sdf_file.c_str()),
			  spec.c_str(),     strlen(spec.c_str()),  
			  &arrays, &e_npart, &e_npart_proc, &e_start);    

	
	    /*
	      if (0 == mpi_rank ) {
	      std::cout << "species: " << spec << " size: " << arrays.l_x << std::endl;
	      std::cout << "npart: " << e_npart << " npart_proc: " << e_npart_proc << " arrays_l: " << arrays.l_x
	      <<   " start: " << e_start << std::endl;
	      
	      for (int k=0; k<arrays.l_x; k++){
	      if ( ( (k%10000) == 0 ) && (mpi_rank == 0) )
	      std::cout << "rank: " << mpi_rank <<  " k: " << k << " x: " << arrays.x[k]
	      << " y: " << arrays.y[k] << " pz: " << arrays.z[k] << std::endl; 
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
	    

	    if (!compression_name.empty()){
	      if (compression_name=="cartesian") {
		
	    // Bining definition for compression
	    int  n_bins[3]  =  {4,4,4};
	    int  p_bins[3]  =  {16,16,16};

            Merger_3d pm(mpi_rank,mpi_size);
	    pm.setVerbose(1);
	    pm.merge(arrays, n_bins, p_bins);
	    
	    // Get mask indexes
	    std::vector<int> vec_mask = pm.get_mask_indexes();
	    if (0 == mpi_rank ) {
	      std::cout << " Merging: rank:"
			<< mpi_rank << " npart: "
			<< arrays.l_px << " tagged indexes: " << vec_mask.size() << std::endl;
	    }

	    // First estimate the updated size (count)
	    int part_size=arrays.l_px;
	    int count{0};	    
	 
	    for (int part_index=0;part_index<part_size; part_index++){
	      bool skip_index=false;
	      if (0 == mpi_rank ) std::cout << " part_index: " <<  part_index << " part_size: " << part_size << std::endl;	    
	      for (size_t mask_index=0;mask_index<vec_mask.size(); mask_index++){
		if (part_index==vec_mask[mask_index]) {skip_index=true; break;}		
	      }//!mask_index
	      if (skip_index==true) continue;
	      //if (0 == mpi_rank ) std::cout << " part_index: " <<  part_index << " accepted " << std::endl;	    
	      count++;
	    }//!part_index

	    if (0 == mpi_rank ) std::cout << "1 count:" <<  count << std::endl;	    

	    // Create buffer
            double array_x[count];
	    double array_y[count];
	    double array_z[count];	    
            double array_w[count];
	    double array_px[count];
            double array_py[count];
	    double array_pz[count];

	    if (0 == mpi_rank ) std::cout << "2" << std::endl;	    	    
	    count=0;
	    for (int part_index=0;part_index<part_size; part_index++){
	      bool skip_index=false;
	      for (size_t mask_index=0;mask_index<vec_mask.size(); mask_index++){
		if (part_index==vec_mask[mask_index]) {skip_index=true; break;}		
	      }//!mask_index
	      if (skip_index==true) continue;
	      array_x[count]=arrays.x[part_index];
	      array_y[count]=arrays.y[part_index];
	      array_z[count]=arrays.y[part_index];	      
	      array_w[count]=arrays.w[part_index];
	      array_px[count]=arrays.px[part_index];
	      array_py[count]=arrays.py[part_index];
	      array_pz[count]=arrays.pz[part_index];	            
	      count++;
	    }//!part_index

	    if (0 == mpi_rank ) std::cout << "3 count: " << count << std::endl;	    	    

	    
	    // Get MPI know the reduction
	    int ntracks_proc=count;
	    int ntracks[mpi_size];	    
	    MPI_Barrier(MPI_COMM_WORLD);	    	    	    
	    MPI_Allgather(&ntracks_proc, 1, MPI_INT,  ntracks, 1, MPI_INT,  MPI_COMM_WORLD);

	    if ( mpi_rank == 0 ){
	      std::cout << std::endl;  
	      for (int i=0 ; i<mpi_size ; i++) {std::cout << " rank: " << i << " ntracks: " << ntracks[i];}  
	      std::cout << std::endl;  
	    }

	    // New total particles  ?
            e_npart=0;
	    for (int i=0; i<mpi_size; i++) e_npart+=ntracks[i];	    
	    
	    // Create Particle species
	    ParticleSpecies e = series.iterations[sdf_iter].particles[spec.c_str()];	    
	    // Create Dataset
	    Datatype datatype = determineDatatype<double>();
	    Extent global_extent = {e_npart};
	    Dataset dataset = Dataset(datatype, global_extent);
	    
	    if (0 == mpi_rank)
	      cout << "Prepared a Dataset of size " << dataset.extent[0] 
		   << " and Datatype " << dataset.dtype
		   << '\n';    

	    // Recalculate particle distribution/proc
            e_start=0;
	    for (int i=0; i<mpi_rank; i++) e_start=e_start+ntracks[i];
	      
	    Offset chunk_offset = {e_start};
	    Extent chunk_extent = {ntracks[mpi_rank]};
           	    
	    // Add reduced particules	    
	    // P
	    if (arrays.l_px>0){ 
	      e["momentum"]["x"].resetDataset(dataset);
	      e["momentum"]["x"].storeChunkRaw( (array_px), chunk_offset, chunk_extent);
	    }
	    
	    if (arrays.l_py>0){ 
	      e["momentum"]["y"].resetDataset(dataset);
	      e["momentum"]["y"].storeChunkRaw( (array_py), chunk_offset, chunk_extent);
	    }
	    
	    if (arrays.l_pz>0){     
	      e["momentum"]["z"].resetDataset(dataset);
	      e["momentum"]["z"].storeChunkRaw( (array_pz), chunk_offset, chunk_extent);
	    }
	    
	    // X
	    if (arrays.l_x>0){         
	      e["position"]["x"].resetDataset(dataset);
	      e["position"]["x"].storeChunkRaw( (array_x), chunk_offset, chunk_extent);
	    }
	    
	    if (arrays.l_y>0){         
	      e["position"]["y"].resetDataset(dataset);
	      e["position"]["y"].storeChunkRaw( (array_y), chunk_offset, chunk_extent);
	    }

	    if (arrays.l_z>0){         
	      e["position"]["z"].resetDataset(dataset);
	      e["position"]["z"].storeChunkRaw( (array_z), chunk_offset, chunk_extent);
	    }	    
	    
	    auto const scalar = openPMD::RecordComponent::SCALAR;
	    
	    if (arrays.l_w>0){         
	      e["weighting"][scalar].resetDataset(dataset);
	      e["weighting"][scalar].storeChunkRaw ( (array_w), chunk_offset, chunk_extent);
	    }
	    
	    // Sync to Disk
	    series.flush();
	      }//! cartesian compression	      
	    }//! compression requested
	    else{
	      // No Particles compression

	      // Create Particle species
	      ParticleSpecies e = series.iterations[sdf_iter].particles[spec.c_str()];
	    
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
	      
	      // Y
	      if (arrays.l_y>0){         
		e["position"]["y"].resetDataset(dataset);
		e["position"]["y"].storeChunkRaw( (arrays.y), chunk_offset, chunk_extent);
	      }
	      
	      // Z
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
	      
	      // Sync to Disk
	      series.flush();
	      
	    }//!else
	  } // species++
	}// ! sdf_files
      return 0;
    } //! Converter3D
    
  }// ! ns_dim_3
} //! ns_converter

using namespace converter::dim_3;

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
    std::cout <<"sdf2opmd_1d:: MPI initialized with size: " << mpi_size << " : " << mpi_rank << std::endl;
    std::cout << std::endl;
  }

  Converter3D opmd_3d;
  opmd_3d.sdf_io(argc, argv);
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
    


    
