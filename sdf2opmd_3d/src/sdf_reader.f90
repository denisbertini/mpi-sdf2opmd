! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE read_support_3d

  USE mpi
  USE sdf

  IMPLICIT NONE

  PUBLIC :: mpi_setup_for_read, read_field_data_r8, read_particle_data

  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18
  INTEGER, PARAMETER :: num = r8
  INTEGER, PARAMETER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER, PARAMETER :: c_max_string_length = 64
  INTEGER, PARAMETER :: c_ndims = 3

  REAL(num), DIMENSION(:,:), POINTER :: current_array
  INTEGER :: current_var, rank

CONTAINS

  SUBROUTINE mpi_setup_for_read

    ! Do the basic MPI setup ops
    INTEGER :: ierr

    !CALL MPI_Init(ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  END SUBROUTINE mpi_setup_for_read



  SUBROUTINE mpi_finish

    INTEGER :: ierr

    !CALL MPI_Finalize(ierr)

  END SUBROUTINE mpi_finish



  SUBROUTINE dims_for_rank(global_dims, rank, n_procs, local_sizes, &
                           local_starts)

    ! Get dimensions of global array belonging to this rank
    INTEGER, DIMENSION(:), INTENT(IN) :: global_dims
    INTEGER, INTENT(IN) :: rank, n_procs
    INTEGER, DIMENSION(:), INTENT(OUT) :: local_sizes, local_starts
    INTEGER, DIMENSION(:), ALLOCATABLE :: dims
    INTEGER :: ierr, i, nproc, sz

    ! Work out this processors piece of domain
    ALLOCATE(dims(c_ndims))

    dims = 0
    CALL MPI_Dims_create(n_procs, c_ndims, dims, ierr)

    dims(1) = n_procs
    dims(2) = 1
    dims(3) = 1 
    
    local_sizes = 0
    local_starts = 0

    nproc = dims (1)
    sz = global_dims(1)/nproc
    
    ! last processor takes excess cells 
    if (rank == (n_procs - 1)) then
       local_starts(1) = sz * rank
       local_sizes(1) = global_dims(1) - max(sz * (nproc - 1), 0)
    else
       local_sizes(1) = sz
       local_starts(1) = sz * rank
    end if
    
    nproc = dims (2)
    sz = global_dims(2)/nproc
    local_sizes(2) = sz
    local_starts(2) = 0

    nproc = dims (3)
    sz = global_dims(3)/nproc
    local_sizes(3) = sz
    local_starts(3) = 0        
    
    DEALLOCATE(dims)

  END SUBROUTINE



  SUBROUTINE create_field_types(n_global, n_local, starts, mpi_basetype, &
                                mpitype, mpi_noghost)

    ! Create the MPI types for fields, without ghost cells
    INTEGER, DIMENSION(:), INTENT(IN) :: n_global, n_local, starts
    INTEGER, INTENT(IN) :: mpi_basetype
    INTEGER, INTENT(OUT) :: mpitype, mpi_noghost
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_starts
    INTEGER :: errcode

    ALLOCATE(local_starts(c_ndims))

    local_starts = 0

    ! MPI type for processors section of global array
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_global, n_local, &
        starts, MPI_ORDER_FORTRAN, mpi_basetype, mpitype, errcode)
    CALL MPI_Type_commit(mpitype, errcode)

    ! MPI type which would be used if we had ghost cells to add/remove
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_local, n_local, &
        local_starts, MPI_ORDER_FORTRAN, mpi_basetype, mpi_noghost, errcode)
    CALL MPI_Type_commit(mpi_noghost, errcode)

    DEALLOCATE(local_starts)

  END SUBROUTINE




    subroutine read_field_data_r8(filename, block_id, field_data, grid_x, grid_y, grid_z,name, &
      units, dims, stagger,  mesh_id, norm, local_sizes, local_starts)
    
    USE sdf_job_info
    ! Read field data in double precision
    ! I check the datatype, but do nothing to convert it if wrong, instead
    ! returning

    character(len=c_max_string_length), intent(in) :: filename
    character(len=c_id_length), intent(in) :: block_id

    ! swap these to single precision for data in single precision
    real(num), dimension(:,:,:), allocatable, intent(out) :: field_data
    real(num), dimension(:), allocatable, intent(out) :: grid_x, grid_y, grid_z
    real(num),intent(out)                             :: norm
    integer, dimension(4), intent(out) :: dims
    character(len=c_id_length), intent(out) :: units, mesh_id
    character(len=c_max_string_length), intent(out) :: name
    integer, intent(out)  :: stagger
    
    character(len=c_id_length) :: code_name, id
    real(num) :: time
    integer, dimension(4), intent(out)  :: local_sizes, local_starts

    integer :: blocktype, datatype, code_io_version, step
    integer :: mpitype, mpi_noghost
    integer :: nblocks, ndims, string_len, total_procs, ierr, i, j
    logical :: restart_flag, found
    type(sdf_file_handle) :: sdf_handle
    type(jobid_type) :: jobid

    CALL sdf_open(sdf_handle, filename, MPI_COMM_WORLD, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      PRINT*, 'Input file contains', nblocks, 'blocks'
    END IF

    CALL sdf_read_blocklist(sdf_handle)

    found = sdf_find_block_by_id(sdf_handle, block_id)

    IF (.NOT. found) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Block not found: ', TRIM(block_id)
      END IF
      STOP
    END IF

    CALL sdf_read_block_header(sdf_handle, id, name, blocktype, ndims, datatype)

    IF (blocktype /= c_blocktype_plain_variable) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Wrong blocktype'
      END IF
      STOP
    END IF

    ! Single precision data has datatype c_datatype_real4
    IF (datatype /= c_datatype_real8) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Wrong datatype'
      END IF
      STOP
    END IF

    ! Check for the correct dimensionality
    IF (ndims /= c_ndims) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Wrong number of dimensions'
      END IF
      STOP
    END IF


    call sdf_read_plain_variable_info(sdf_handle, dims, units, mesh_id, stagger, norm)    


    if (rank == 0) then
       print*, ' read_field_data_r8():' 
       print*, '               type: ', datatype
       print*, '              ndims: ', ndims
       print*, '               dims: ', dims(1:ndims)
    end if

    CALL MPI_Comm_size(MPI_COMM_WORLD, total_procs, ierr)
    CALL dims_for_rank(dims, rank, total_procs, local_sizes, local_starts)

    print*, ' read_field_data_r8():' 
    print*, 'rank: ', rank, ' sizes: ', local_sizes, ' starts: ', local_starts    

    ALLOCATE(field_data(local_sizes(1), local_sizes(2), local_sizes(3)))
    field_data = 0.0_num

    CALL create_field_types(dims, local_sizes, local_starts, mpireal, &
        mpitype, mpi_noghost)

    CALL sdf_read_plain_variable(sdf_handle, field_data, mpitype, mpi_noghost)

    IF (rank == 0) PRINT*, 'Grid name is ', mesh_id

    ALLOCATE(grid_x(dims(1)+1))
    ALLOCATE(grid_y(dims(2)+1))
    ALLOCATE(grid_z(dims(3)+1))

    grid_x = 0.0_num
    grid_y = 0.0_num
    grid_z = 0.0_num

    found = sdf_find_block_by_id(sdf_handle, mesh_id)

    IF (.NOT. found) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Block not found: ', TRIM(mesh_id)
      END IF
      STOP
    END IF

    CALL sdf_read_srl_plain_mesh(sdf_handle, grid_x, grid_y, grid_z)

    CALL sdf_close(sdf_handle)

    CALL MPI_Type_free(mpitype, ierr)
    CALL MPI_Type_free(mpi_noghost, ierr)

  END SUBROUTINE read_field_data_r8



  SUBROUTINE particles_for_rank(npart, rank, n_procs, npart_proc, start)

    ! Work out this ranks number of processors
    INTEGER(i8), INTENT(IN) :: npart
    INTEGER, INTENT(IN) :: rank, n_procs
    INTEGER(i8), INTENT(OUT) :: npart_proc, start
    INTEGER(i8) :: npart_proc_gen

    ! Work out this processors number of particles
    ! Have to handle non divisibility

    ! Size per processor
    npart_proc_gen = npart / n_procs

    IF (rank == n_procs) THEN
      ! Last processor takes excess particles to make starts easier
      npart_proc = npart - MAX(npart_proc_gen * (n_procs - 1), 0_i8)
    ELSE
      npart_proc = npart_proc_gen
    END IF

    start = rank * npart_proc_gen + 1

  END SUBROUTINE



  SUBROUTINE read_particle_data(filename, species_name, particle_data, npart_global, npart_proc, start, columns)

    USE sdf_job_info
    ! Read all particle data for named species into data array
    ! First three columns will be 3-D grid data
    ! Next will be one per particle variable
    CHARACTER(LEN=c_max_string_length), INTENT(IN) :: filename, species_name
    ! Swap these to single precision for data in single precision
    ! Read however many particle variables are present
    ! Target property so we can properly use the iterators
    REAL(num), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT) :: particle_data
    INTEGER(i8)                                   , INTENT(OUT) :: npart_global, npart_proc, start
    INTEGER(i8)                                                 :: npart 
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, found_mesh_id, species_id
    CHARACTER(LEN=c_id_length), dimension(:), allocatable, target, intent(out) :: columns
    CHARACTER(LEN=c_id_length) :: units
    CHARACTER(LEN=c_max_string_length) :: name
    REAL(num) :: time
    INTEGER :: nblocks, ndims, string_len, total_procs, ierr
    INTEGER :: vars_per_species, iblock, step, geometry
    INTEGER :: blocktype, datatype, code_io_version, mpitype
    LOGICAL :: restart_flag, found
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(jobid_type) :: jobid

    CALL sdf_open(sdf_handle, filename, MPI_COMM_WORLD, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      PRINT*, 'Input file contains', nblocks, 'blocks'
    END IF

    CALL sdf_read_blocklist(sdf_handle)

    vars_per_species = 0

    ! Run through file once to see how many particle variables there are
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype /= c_blocktype_point_variable) CYCLE

      CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
          units, species_id)

      ! To read only selected variables, select the names here
      IF (TRIM(species_id) == species_name) THEN
         vars_per_species = vars_per_species + 1
         found_mesh_id = mesh_id
         npart_global = npart         
      END IF
    END DO

    IF (vars_per_species < 1) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'No variables found for species: ', TRIM(species_name)
      END IF
      STOP
    END IF

    ! Get number of particles
    found = sdf_find_block_by_id(sdf_handle, found_mesh_id)

    IF (.NOT. found) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Block not found: ', TRIM(found_mesh_id)
      END IF
      STOP
    END IF

    CALL sdf_read_point_mesh_info(sdf_handle, npart)

    IF (npart <= 0) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'No particles for species: ', TRIM(species_name)
      END IF
      STOP
    END IF

   IF (rank == 0) PRINT*, 'Found  species_name: ', species_name, ' npart: ', npart, &
        ' with: ',  vars_per_species, ' particle variables'    

    CALL MPI_Comm_size(MPI_COMM_WORLD, total_procs, ierr)
    CALL particles_for_rank(npart, rank, total_procs, npart_proc, start)

    ! Allocate arrays
    ALLOCATE(particle_data(npart_proc, vars_per_species+c_ndims))

    CALL sdf_seek_start(sdf_handle)
    CALL create_particle_type(npart_proc, total_procs, mpitype)

    ! Setup globals for access within it_part
    current_array => particle_data

    ! Read the grid
    current_var = 1
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype /= c_blocktype_point_mesh) CYCLE

      CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)

      IF (TRIM(species_id) == species_name) THEN
        ! Read this procs grid
        CALL sdf_read_point_mesh(sdf_handle, npart_proc, mpitype, it_part_mesh)
        IF (rank == 0) PRINT*, 'Columns 1 to',c_ndims,' Columns', current_var, ' is ', block_id
      END IF
    END DO

    ! Rewind and read the data
    ! This is inefficient but simplest
    CALL sdf_seek_start(sdf_handle)
    current_var = c_ndims + 1
    ! Read the data
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype /= c_blocktype_point_variable) CYCLE

      CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
          units, species_id)

      ! To read only selected variables, select the names here.
      ! These must match line 305
      IF (TRIM(species_id) == species_name) THEN
        ! Read this procs particles
        CALL sdf_read_point_variable(sdf_handle, npart_proc, mpitype, &
            it_part_arr)
        IF (rank == 0) PRINT*, 'Column ', current_var, ' is ', block_id
        current_var = current_var + 1
      END IF
    END DO

    ! Additional particle info
    allocate (columns(current_var));
    columns(1) = 'gridx'
    columns(2) = 'gridy'
    columns(3) = 'gridz'    
       
    CALL sdf_seek_start(sdf_handle)
    current_var = c_ndims+1

    ! Read again the data
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype /= c_blocktype_point_variable) CYCLE

      CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
          units, species_id)

      IF (TRIM(species_id) == species_name) THEN
        !IF (rank == 0) PRINT*, 'AGAIN Column ', current_var, ' is ', block_id
        columns(current_var) = block_id
        current_var = current_var + 1        
      END IF
    END DO
    
    CALL sdf_close(sdf_handle)
    CALL MPI_Type_free(mpitype, ierr)

  END SUBROUTINE read_particle_data



  ! Particle iterator borrowed from setup and adapted

  FUNCTION it_part_mesh(array, npart_this_it, start, direction, param)

    ! Iterator for reading grid
    ! Direction denotes x,y,z
    REAL(num) :: it_part_mesh
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER(i8) :: ipart
    INTEGER(i8), SAVE :: ipart_done

    IF (start) ipart_done = 0

    DO ipart = 1, npart_this_it
      current_array(ipart+ipart_done, current_var+direction-1) = array(ipart)
    END DO

    ipart_done = ipart_done + npart_this_it
    it_part_mesh = 0

  END FUNCTION it_part_mesh



  ! Particle iterator borrowed from setup and adapted

  FUNCTION it_part_arr(array, npart_this_it, start, param)

    ! Iterator for particle variable
    REAL(num) :: it_part_arr
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER(i8) :: ipart
    INTEGER(i8), SAVE :: ipart_done

    IF (start) ipart_done = 0

    DO ipart = 1, npart_this_it
      current_array(ipart+ipart_done, current_var) = array(ipart)
    END DO

    ipart_done = ipart_done + npart_this_it
    it_part_arr = 0

  END FUNCTION it_part_arr



  ! Adapted from create_particle_subtypes in mpi_subtype_control

  SUBROUTINE create_particle_type(npart_in, nproc, subtype)

    ! Type for particle reading. Must account for where in data to start
    INTEGER(i8), INTENT(IN) :: npart_in
    INTEGER, INTENT(IN) :: nproc
    INTEGER, INTENT(OUT) :: subtype
    INTEGER(i8), DIMENSION(1) :: npart_local
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_each_rank
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(3) :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: particles_to_skip, total_particles
    INTEGER :: i, mpitype, basetype, typesize, errcode

    npart_local = npart_in

    ALLOCATE(npart_each_rank(nproc))

    ! Create the subarray for the particles in this problem: subtype decribes
    ! where this process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, &
        npart_each_rank, 1, MPI_INTEGER8, MPI_COMM_WORLD, errcode)

    particles_to_skip = 0
    DO i = 1, rank
      particles_to_skip = particles_to_skip + npart_each_rank(i)
    END DO

    total_particles = particles_to_skip
    DO i = rank + 1, nproc
      total_particles = total_particles + npart_each_rank(i)
    END DO

    DEALLOCATE(npart_each_rank)

    basetype = mpireal
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)

    ! If npart_in is bigger than an integer then the data will not
    ! get written/read properly. This would require about 48GB per processor
    ! so it is unlikely to be a problem any time soon.
    lengths(1) = 1
    lengths(2) = INT(npart_in)
    lengths(3) = 1
    disp(1) = 0
    disp(2) = particles_to_skip * typesize
    disp(3) = total_particles * typesize
    types(1) = MPI_LB
    types(2) = basetype
    types(3) = MPI_UB

    mpitype = 0
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    subtype = mpitype

  END SUBROUTINE create_particle_type



  subroutine read_derived_variables(filename)
    
    use sdf_job_info
    
    character(len=c_max_string_length), intent(in) :: filename
    character(len=c_id_length) :: block_id

    integer, dimension(4) :: dims
    character(len=c_max_string_length) :: name
    integer :: stagger
    
    character(len=c_id_length) :: code_name, id
    real(num) :: time
    integer, dimension(4) :: local_sizes, local_starts

    integer :: blocktype, datatype, code_io_version, step
    integer :: mpitype, mpi_noghost
    integer :: nblocks, ndims, string_len, total_procs, ierr, i, j
    logical :: restart_flag, found
    integer :: iblock
    character(len=c_id_length) :: tmpblock_id

    type(sdf_file_handle) :: sdf_handle
    type(jobid_type) :: jobid

    call sdf_open(sdf_handle, filename, mpi_comm_world, c_sdf_read)   
    call sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
         string_len, restart_flag)
    
    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)
    
    call sdf_read_blocklist(sdf_handle)
  
    
    if (rank == 0) print*, "Input file contains: ", nblocks, " blocks"

    tmpblock_id=" "
    do iblock = 1, nblocks
       call sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
            ndims, datatype)
       if (rank == 0 ) print*,'Block info: ', iblock, " : ",  name, " : ",  block_id, " : " , blocktype
       
      ! selection       
       if (index(name, 'Derived') > 0) then  
          tmpblock_id = trim(tmpblock_id) // ':' // trim(block_id)
       end if
    end do

    
  end subroutine read_derived_variables


  
END MODULE read_support_3d

