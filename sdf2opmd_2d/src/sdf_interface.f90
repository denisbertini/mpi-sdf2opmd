module sdf_interface
use, intrinsic :: iso_c_binding

use read_support
use sdf

implicit none

contains
 
  subroutine init_read() bind(C)
    call mpi_setup_for_read
  end subroutine init_read
  
  subroutine read_particle(cstring, clen, cspecies, clenspec, arrays, npart, npart_proc, start) bind(C)
    integer(c_int), intent(in), value :: clen
    character(c_char), intent(in) :: cstring(clen)    
    integer(c_int), intent(in), value :: clenspec
    character(c_char), intent(in) :: cspecies(clenspec)    
    character(len=clen) :: fstring
    character(len=clenspec) :: fspecies    
    integer(c_long), intent(out) :: npart, npart_proc, start 

    
    real(c_double), dimension(:,:), allocatable :: field_data
    real(c_double), dimension(:,:), allocatable, target  :: particle_data
    real(c_double), dimension(:),   allocatable :: grid_x, grid_y
    real(c_double), allocatable, target, save          :: eta_px(:), eta_py(:), eta_pz(:)  
    character(LEN=c_max_string_length) :: read_dir = './Data/'
    character(LEN=c_max_string_length) :: filename, species_name
    character(LEN=c_id_length) :: block_id
    integer                    :: i , j
    
    type, bind(C) :: part
       !integer (c_int) :: lenc
       !type(c_ptr)     :: c        
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
