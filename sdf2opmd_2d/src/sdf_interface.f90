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
    character(LEN=c_id_length),   dimension(:),   allocatable, target  :: columns
    real(c_double), dimension(:),   allocatable          :: grid_x, grid_y
    real(c_double), allocatable, target, save            :: eta_px(:), eta_py(:), eta_pz(:)
    real(c_double), allocatable, target, save            :: eta_x(:), eta_y(:), eta_w(:)
    real(c_double), allocatable, target, save            :: eta_vx(:), eta_vy(:), eta_vz(:)
    real(c_double), allocatable, target, save            :: eta_ek(:), eta_rm(:), eta_gm(:)

    character(LEN=c_max_string_length) :: read_dir = './Data/'
    character(LEN=c_max_string_length) :: filename, species_name
    character(LEN=c_id_length) :: block_id
    integer                    :: i , j, ipos
    
    type, bind(C) :: part
       integer (c_int) :: l_x, l_y, l_w, l_px, l_py, l_pz, l_vx, l_vy, l_vz, l_ek, l_rm, l_gm
       type(c_ptr)     ::   x,   y,   w,   px,   py,   pz,   vx,   vy,   vz,   ek,   rm,   gm 
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
    call read_particle_data(filename, species_name, particle_data, npart, npart_proc, start, columns)
    ! Example calc - calculate average column 3 per processor
    !print*, 'deuteron average px on rank ', rank, ' is ', &
    !     sum(particle_data(:,4)) / size(particle_data(:,4))
    !print*, 'columns info on rank: ', rank, ' size: ', &
    !     size(columns), columns

    
    ! Convert particles informations according to the origin SDF input
    if ( get_index(columns, 'gridx', species_name) > 0) then 
       arrays%l_x  = size(particle_data(:, get_index(columns, 'gridx', species_name))) 
       if (allocated(eta_x)) deallocate(eta_x)  
       allocate (eta_x(arrays%l_x))
       eta_x(1:arrays%l_x) = particle_data(1:arrays%l_x, get_index(columns, 'gridx', species_name)) 
       arrays%x = c_loc(eta_x)
    else
       arrays%l_x = 0
       arrays%x   = c_null_ptr         
    end if

    if ( get_index(columns, 'gridy', species_name) > 0) then 
       arrays%l_y  = size(particle_data(:, get_index(columns, 'gridy', species_name))) 
       if (allocated(eta_y)) deallocate(eta_y)  
       allocate (eta_y(arrays%l_y))
       eta_y(1:arrays%l_y) = particle_data(1:arrays%l_y, get_index(columns, 'gridy', species_name)) 
       arrays%y = c_loc(eta_y)
    else
       arrays%l_y = 0
       arrays%y   = c_null_ptr         
    end if

    if ( get_index(columns, 'weight', species_name) > 0) then
       arrays%l_w  = size(particle_data(:, get_index(columns, 'weight', species_name)))
       if (allocated(eta_w)) deallocate(eta_w)  
       allocate (eta_w(arrays%l_w))       
       eta_w(1:arrays%l_w) = particle_data(1:arrays%l_w, get_index(columns, 'weight', species_name))
       arrays%w = c_loc(eta_w)       
    else
       arrays%l_w = 0
       arrays%w   = c_null_ptr         
    end if

    if ( get_index(columns, 'px', species_name) > 0) then
       arrays%l_px = size(particle_data(:, get_index(columns, 'px', species_name)))
       if (allocated(eta_px)) deallocate(eta_px)  
       allocate (eta_px(arrays%l_px))
       eta_px(1:arrays%l_px) = particle_data(1:arrays%l_px, get_index(columns, 'px', species_name))
       arrays%px = c_loc(eta_px)       
    else
       arrays%l_px = 0
       arrays%px   = c_null_ptr         
    end if
    
    if ( get_index(columns, 'py', species_name) > 0) then
       arrays%l_py = size(particle_data(:, get_index(columns, 'py', species_name)))
       if (allocated(eta_py)) deallocate(eta_py)  
       allocate (eta_py(arrays%l_py))       
       eta_py(1:arrays%l_py) = particle_data(1:arrays%l_py, get_index(columns, 'py', species_name))
       arrays%py = c_loc(eta_py)       
    else
       arrays%l_py = 0
       arrays%py   = c_null_ptr         
    end if

    if ( get_index(columns, 'pz', species_name) > 0) then
       arrays%l_pz = size(particle_data(:, get_index(columns, 'pz', species_name)))
       if (allocated(eta_pz)) deallocate(eta_pz)  
       allocate (eta_pz(arrays%l_pz))
       eta_pz(1:arrays%l_pz) = particle_data(1:arrays%l_pz, get_index(columns, 'pz', species_name)) 
       arrays%pz = c_loc(eta_pz)
    else
       arrays%l_pz = 0
       arrays%pz   = c_null_ptr         
    end if

    if ( get_index(columns, 'vx', species_name) > 0) then
       arrays%l_vx = size(particle_data(:, get_index(columns, 'vx', species_name)))
       if (allocated(eta_vx)) deallocate(eta_vx)  
       allocate (eta_vx(arrays%l_vx))       
       eta_vx(1:arrays%l_vx) = particle_data(1:arrays%l_vx, get_index(columns, 'vx', species_name)) 
       arrays%vx = c_loc(eta_vx)
    else
       arrays%l_vx = 0
       arrays%vx   = c_null_ptr         
    end if
    
    if ( get_index(columns, 'vy', species_name) > 0) then
       arrays%l_vy = size(particle_data(:, get_index(columns, 'vy', species_name))) 
       if (allocated(eta_vy)) deallocate(eta_vy)  
       allocate (eta_vy(arrays%l_vy))
       eta_vy(1:arrays%l_vy) = particle_data(1:arrays%l_vy, get_index(columns, 'vy', species_name)) 
       arrays%vy = c_loc(eta_vy)
    else
       arrays%l_vy = 0
       arrays%vy   = c_null_ptr         
    end if
    
    if ( get_index(columns, 'vz', species_name) > 0) then
       arrays%l_vz = size(particle_data(:, get_index(columns, 'vz', species_name)))        
       if (allocated(eta_vz)) deallocate(eta_vz)  
       allocate (eta_vz(arrays%l_vz))
       eta_vz(1:arrays%l_vz) = particle_data(1:arrays%l_vz, get_index(columns, 'vz', species_name)) 
       arrays%vz = c_loc(eta_vz)
    else
       arrays%l_vz = 0
       arrays%vz   = c_null_ptr         
    end if
    
    if ( get_index(columns, 'ek', species_name) > 0) then
       arrays%l_ek = size(particle_data(:, get_index(columns, 'ek', species_name))) 
       if (allocated(eta_ek)) deallocate(eta_ek)  
       allocate (eta_ek(arrays%l_ek))
       eta_ek(1:arrays%l_ek) = particle_data(1:arrays%l_ek, get_index(columns, 'ek', species_name)) 
       arrays%ek = c_loc(eta_ek)
    else
       arrays%l_ek = 0
       arrays%ek   = c_null_ptr                
    end if
    
    if ( get_index(columns, 'relativistic mass', species_name) > 0) then
       arrays%l_rm = size(particle_data(:, get_index(columns, 'relativistic mass', species_name))) 
       if (allocated(eta_rm)) deallocate(eta_rm)  
       allocate (eta_rm(arrays%l_rm))
       eta_rm(1:arrays%l_rm) = particle_data(1:arrays%l_rm, get_index(columns, 'relativistic mass', species_name)) 
       arrays%rm = c_loc(eta_rm)
    else
       arrays%l_rm = 0
       arrays%rm   = c_null_ptr                
    end if
    
    if ( get_index(columns, 'gamma', species_name) > 0) then
       arrays%l_gm = size(particle_data(:, get_index(columns, 'gamma', species_name))) 
       if (allocated(eta_gm)) deallocate(eta_gm)  
       allocate (eta_gm(arrays%l_gm))
       eta_gm(1:arrays%l_vz) = particle_data(1:arrays%l_gm, get_index(columns, 'gamma', species_name)) 
       arrays%gm = c_loc(eta_gm)        
    else
       arrays%l_gm = 0
       arrays%gm   = c_null_ptr                
    end if
    
    ! Release memory ONLY from fortran side 
    deallocate(particle_data, columns)

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

  function get_index(cols, token, spec_name) result(ipos)
    character(len=c_id_length), dimension(:), intent(in)  :: cols
    character(len=*), intent(in)                          :: token
    character(len=c_max_string_length), intent(in)        :: spec_name    
    integer                                               :: i, ipos
    
    do i = 1, size(cols)
       if ((index(cols(i), trim(token) )) > 0 ) then
          ipos = i
          !print*, token , ' found at column', ipos
          return
       end if   
    end do

    print*, 'Warning variable token:', trim(token), ' not matched for species: ', spec_name    
    ipos = -1
    
  end function get_index
  
end module sdf_interface
