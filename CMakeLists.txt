cmake_minimum_required (VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_FLAGS "-Wall -O3 -std=c11 ${CMAKE_C_FLAGS}")
project (JellyBelly)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}"
                    "klib"
                    "spaces")

# add the executable
add_executable(JellyBelly src/main.c
                          src/bellyFun.c
                          src/bellyHash.c
                          src/bellyMisc.c)

add_executable(create_spaces spaces/create_spaces.c
                             spaces/space_functions.c)

add_executable(distMatrix utilities/distMatrix.c)

find_library(ZLIB_LIBRARY
             NAMES z)
find_library(MLIB_LIBRARY
             NAMES m)


get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
target_link_libraries(JellyBelly
                      ${ZLIB_LIBRARY}
                      ${MLIB_LIBRARY})

target_link_libraries(create_spaces
                      ${MLIB_LIBRARY})

target_link_libraries(distMatrix
                      ${MLIB_LIBRARY})
