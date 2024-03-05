# -- CMake configuration -------------------------------------------------------
cmake_minimum_required(VERSION 3.24)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# -- Project description -------------------------------------------------------
project(LRAP 
    LANGUAGES Fortran)

# -- Languages -----------------------------------------------------------------
enable_language(Fortran)

# -- Set up build and additional directories -----------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# -- Subprojects ---------------------------------------------------------------
add_subdirectory(MARIA-Fortran)

# -- Source files --------------------------------------------------------------
add_subdirectory(src)

# -- Set up unit tests ---------------------------------------------------------
include(CTest)
enable_testing()
