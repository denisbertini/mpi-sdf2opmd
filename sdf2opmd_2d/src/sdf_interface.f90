module sdf_interface
  use, intrinsic :: iso_c_binding
  use read_support
  use sdf
  
  implicit none

contains
 
  subroutine init_read() bind(C)
    call mpi_setup_for_read
  end subroutine init_read

  subroutine read_fields(cstring, clen, blockid , blen, cfield_data_x) bind(C)
    
    integer(c_int), intent(in), value                :: clen
    character(c_char), intent(in)                    :: cstring(clen)
    integer(c_int), intent(in), value                :: blen
    character(c_char), intent(in)                    :: blockid(blen)    
    character(len=clen)                              :: fstring
    character(len=blen)                              :: cblockid

    real(c_double), dimension(:,:), allocatable      :: field_data
    real(c_double)                                   :: stag_x, stag_y
    character(LEN=c_id_length)                       :: block_id
    real(num), dimension(:),    allocatable          :: grid_x, grid_y
    real(num)                                        :: norm
    character(len=c_id_length)                       :: units, mesh_id
    character(len=c_max_string_length)               :: name
    integer, dimension(4)                            :: dims, local_starts, local_sizes
    integer                                          :: i, j, stagger_x, stagger_y, l_index 
    character(LEN=c_max_string_length)               :: filename
    real(c_double), allocatable, target, save       :: lin_field_data(:), cgrid_x(:), cgrid_y(:)    

   type, bind(C) :: field
       integer (c_int) :: global_sx, global_sy, l_sx, l_sy, l_dx, l_dy, l_gridx, l_gridy, l_data_size
       real(c_double)  :: staggger
       type(c_ptr)     :: gridx, gridy, l_field_data 
    end type  field

    type(field), intent(out) :: cfield_data_x 

    ! Convert filename 
    fstring  = s_copy(cstring, clen)
    cblockid = s_copy(blockid, blen)
    filename = fstring
    block_id = cblockid  
    
    ! Get Fields Data
    print*, ' block_id: ', block_id
      call read_field_data_r8(filename, block_id, field_data, grid_x, grid_y, name, units, &
           dims, stagger_x, mesh_id, norm, local_sizes, local_starts)

    if (stagger_x == c_stagger_face_x) then
         cfield_data_x%staggger=0.5
    end if

    ! Domain partition 
    cfield_data_x%global_sx = dims(1)
    cfield_data_x%global_sy = dims(2)    
    cfield_data_x%l_sx = local_sizes(1)
    cfield_data_x%l_sy = local_sizes(2)
    cfield_data_x%l_dx = local_starts(1)
    cfield_data_x%l_dy = local_starts(2)

    ! Arrays sizes
    cfield_data_x%l_gridx = size(grid_x)
    cfield_data_x%l_gridy = size(grid_y)
    cfield_data_x%l_data_size = size(field_data)
    
    ! Assign arrays
    if (allocated(cgrid_x)) deallocate(cgrid_x)  
    allocate (cgrid_x(cfield_data_x%l_gridx))
    if (allocated(cgrid_y)) deallocate(cgrid_y)  
    allocate (cgrid_y(cfield_data_x%l_gridy))
    if (allocated(lin_field_data)) deallocate(lin_field_data)  
    allocate (lin_field_data(0:cfield_data_x%l_data_size-1))
    
    cgrid_x(1:cfield_data_x%l_gridx) = grid_x(1:cfield_data_x%l_gridx)
    cfield_data_x%gridx = c_loc(cgrid_x)
    cgrid_y(1:cfield_data_x%l_gridy) = grid_y(1:cfield_data_x%l_gridy)
    cfield_data_x%gridy = c_loc(cgrid_y)

   
    ! Transform field to linear array
     do i=1, local_sizes(1)
       do j=1, local_sizes(2)
          ! offset = i_col*N_rows+i_row 
          l_index = ((j-1) * local_sizes(1)) + (i-1)  
          lin_field_data(l_index) = field_data(i,j)
          !if ( (i<10) .and. (j<2) ) then  
          !   print *,' i: ', i, ' j: ', j, ' lin: ' , l_index, ' lin_data: ', lin_field_data(l_index), &
          !        ' field_data: ', field_data(i,j)  
          !end if
       end do
    end do 

   cfield_data_x%l_field_data = c_loc(lin_field_data)
      
  end subroutine read_fields
    
  subroutine read_particle(cstring, clen, cspecies, clenspec, arrays, npart, npart_proc, start) bind(C)
    integer(c_int), intent(in), value :: clen
    character(c_char), intent(in) :: cstring(clen)    
    integer(c_int), intent(in), value :: clenspec
    character(c_char), intent(in) :: cspecies(clenspec)    
    character(len=clen) :: fstring
    character(len=clenspec) :: fspecies    
    integer(c_long), intent(out) :: npart, npart_proc, start 

    real(c_double), dimension(:,:), allocatable, target  :: particle_data
    real(c_double), dimension(:),   allocatable :: grid_x, grid_y
    real(c_double), allocatable, target, save          :: eta_px(:), eta_py(:), eta_pz(:)  
    character(LEN=c_max_string_length) :: read_dir = './Data/'
    character(LEN=c_max_string_length) :: filename, species_name
    character(LEN=c_id_length) :: block_id
    integer                    :: i , j
    
    type, bind(C) :: part
       integer (c_int) :: l_px, l_py, l_pz
       type(c_ptr)     :: px, py, pz 
    end type  part
    
    type(part), intent(out) :: arrays 
    real (c_double), pointer :: c_array(:) 

    fstring  = s_copy(cstring, clen)
    fspecies = s_copy(cspecies, clenspec)
    filename = fstring
    species_name = fspecies  

    if ( rank == 0 ) then
       print*, 'SDF input filename: ', filename
       print*, 'Fetching data for species: ', species_name
    end if
    
    ! Particle data
    call read_particle_data(filename, species_name, particle_data, npart, npart_proc, start)
    ! Example calc - calculate average column 3 per processor
    !print*, 'deuteron average px on rank ', rank, ' is ', &
    !     sum(particle_data(:,4)) / size(particle_data(:,4))

    ! Get part. array size
    arrays%l_px = size(particle_data(:,4)) 
    arrays%l_py = size(particle_data(:,5))
    arrays%l_pz = size(particle_data(:,6))
    
    ! Allocate an array and make it available in C
    !print*, 'species_name size px, py pz:', arrays%l_px, arrays%l_py, arrays%l_pz
    ! If already static array ETA allocated deallocate  

    if (allocated(eta_px)) deallocate(eta_px)  
    allocate (eta_px(arrays%l_px))
    if (allocated(eta_py)) deallocate(eta_py)  
    allocate (eta_py(arrays%l_py))
    if (allocated(eta_pz)) deallocate(eta_pz)  
    allocate (eta_pz(arrays%l_pz))

    ! Copy contents of particle_data to Eta array
    eta_px(1:arrays%l_px) = particle_data(1:arrays%l_px,4)
    eta_py(1:arrays%l_py) = particle_data(1:arrays%l_py,5)
    eta_pz(1:arrays%l_pz) = particle_data(1:arrays%l_pz,6)
    arrays%px = c_loc(eta_px)
    arrays%py = c_loc(eta_py)
    arrays%pz = c_loc(eta_pz)    
    
    ! Release memory ONLY from fortran side 
    deallocate(particle_data)
  end subroutine read_particle


  
  function s_copy(cstring, clen) result(fstring)   
    integer(c_int), intent(in), value :: clen
    character(c_char), intent(in)     :: cstring(clen)
    character(len=clen)               :: fstring
    integer                           :: j
    
    fstring = ' ' 
    do j = 1, clen
       if (cstring(j) == c_null_char) exit
       fstring(j:j) = cstring(j)
    end do
  end function s_copy
 
end module sdf_interface
