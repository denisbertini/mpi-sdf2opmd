
find_package(openPMD 0.14.5 CONFIG)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wno-narrowing")

SET ( EXECUTABLENAME sdf2opmd_2d )

SET ( ${EXECUTABLENAME}_SRCS
        sdf2opmd_2d.cxx
)

include_directories(
        .
        ${CMAKE_SOURCE_DIR}/sdf2opmd_2d/src
        ${CMAKE_SOURCE_DIR}/sdf2opmd_2d/merger	
)

ADD_EXECUTABLE(${EXECUTABLENAME}
        ${${EXECUTABLENAME}_SRCS}
)

TARGET_LINK_LIBRARIES (${EXECUTABLENAME} PRIVATE openPMD::openPMD SDFOPMD2D "stdc++fs")

ADD_TEST(${EXECUTABLENAME} ${EXECUTABLENAME})

INSTALL ( TARGETS ${EXECUTABLENAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

INCLUDE_DIRECTORIES ()

