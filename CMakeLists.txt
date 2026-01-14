cmake_minimum_required (VERSION 3.8)
project (ParticleGeneratorGPU LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(ROOT CONFIG REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Foam)
add_definitions("-DSMASH_TOP_LEVEL_DIR=\"$ENV{SMASH_DIR}\"")
add_definitions(-DSMASH_USE_ROOT)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
#Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
        message(FATAL_ERROR "You are trying to configure in the source directory!")
endif()

#include_directories(include)
file(GLOB SOURCES "src/*.cxx")

# Add source to this project's executable.
add_executable(ThermalBlastMC ${SOURCES})

find_package(SMASH)
if(SMASH_FOUND)
    message(STATUS "SMASH libraries = ${SMASH_LIBRARIES}")
    include_directories(${SMASH_INCLUDE_DIR})
    add_definitions("-DSMASH_INPUT_DIR=\"${SMASH_INPUT_FILES_DIR}\"")
else()
    message(FATAL_ERROR "SMASH libraries not found!")
endif()
message(STATUS "ROOT libraries = ${ROOT_LIBRARIES}")
target_link_libraries(ThermalBlastMC PUBLIC ${ROOT_LIBRARIES} ${SMASH_LIBRARIES})
target_include_directories(ThermalBlastMC PRIVATE
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/include
)
