


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
	${EPOCH_PATH}/epoch2d/obj
	}
)

link_directories(${OTHER_LINK_DIR} ${SDF_PATH}/lib ${EPOCH_PATH}/epoch2d/obj/)

set( SDFOPMD2D_HEADERS
	src/cxxopts.hpp
)

set( SDFOPMD2D_SRCS
	src/sdf_interface.f90
	src/sdf_reader.f90
)


add_library( SDFOPMD2D SHARED ${SDFOPMD2D_SRCS})

target_link_libraries(SDFOPMD2D
   PRIVATE
   MPI::MPI_C
   ${MPI_Fortran_LIBRARIES}
   ${MPI_LIBRARIES}
   ${SDF_LIB}
   ${CMAKE_SOURCE_DIR}/sdf2opmd_2d/epoch_libs/epoch2d_lib.a	   
   )

set_target_properties(SDFOPMD2D PROPERTIES LINKER_LANGUAGE CXX)

# Add executables
add_subdirectory(converter)
add_subdirectory(reader)

install(
        TARGETS SDFOPMD2D
        LIBRARY DESTINATION lib
)


