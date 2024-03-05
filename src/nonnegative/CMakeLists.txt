add_subdirectory(matrix)
add_subdirectory(tensor)
add_subdirectory(matrix_completion)

# Test with 2D Smolukhovski equation
set_source_files_properties(matrix_smolukh.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_smolukh
    matrix_smolukh.f90)
target_link_libraries(matrix_smolukh
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} 
    maria_la_core nonnegative_matrix)

# Test with QTT mixture of Gaussians
set_source_files_properties(tensor_gaussian.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(tensor_gaussian
    tensor_gaussian.f90)
target_link_libraries(tensor_gaussian
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} 
    maria_la_core nonnegative_tensor)

# Test with matrix completion
set_source_files_properties(matrix_cmplt.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_cmplt
    matrix_cmplt.f90)
target_link_libraries(matrix_cmplt
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} 
    maria_la_core maria_la_utils nn_matrix_completion)

