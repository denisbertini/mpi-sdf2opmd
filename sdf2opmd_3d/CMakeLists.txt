
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "ROOT_INCLUDE_PATH : ${ROOT_INCLUDE_DIR}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ compiler Flags are : ${CMAKE_CXX_FLAGS}")
message(STATUS "C++ compiler identification is : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "EPOCH     PATH is : ${EPOCH_PATH}")
message(STATUS "EPOCH SDF PATH is : ${SDF_PATH}")



include_directories(include
        ${MPI_INCLUDE_PATH}  
	${SDF_PATH}/include
	${EPOCH_PATH}/epoch3d/obj
	src/
	}
)

link_directories(${OTHER_LINK_DIR} ${SDF_PATH}/lib ${EPOCH_PATH}/epoch3d/obj/)

set( SDFOPMD3D_HEADERS
	src/cxxopts.hpp
	src/utils.h
)

set( SDFOPMD3D_SRCS
	src/sdf_interface.f90
	src/sdf_reader.f90
	merger/merger_3d.cxx
)


add_library( SDFOPMD3D SHARED ${SDFOPMD3D_SRCS})

target_link_libraries(SDFOPMD3D
   PRIVATE
   MPI::MPI_C
   ${MPI_Fortran_LIBRARIES}
   ${MPI_LIBRARIES}
   ${SDF_LIB}
   ${CMAKE_SOURCE_DIR}/sdf2opmd_3d/epoch_libs/epoch3d_lib.a	   
   )

set_target_properties(SDFOPMD3D PROPERTIES LINKER_LANGUAGE CXX)

#apply_global_cxx_flags_to_all_targets()
#remove_flag_from_target(SDFOPMD3D -Wmissing-include-dirs)

# Add executables
add_subdirectory(converter)

install(
        TARGETS SDFOPMD3D
        LIBRARY DESTINATION lib
)


