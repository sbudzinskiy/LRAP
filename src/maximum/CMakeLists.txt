add_subdirectory(matrix)
add_subdirectory(tensor)

# Test with random orthogonal matrix
set_source_files_properties(matrix_randort.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_randort
    matrix_randort.f90)
target_link_libraries(matrix_randort
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_matrix)

# Test with low-rank random orthogonal matrix
set_source_files_properties(matrix_lr_randort.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_lr_randort
    matrix_lr_randort.f90)
target_link_libraries(matrix_lr_randort
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_matrix)

# Test with the identity matrix
set_source_files_properties(matrix_identity.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_identity
    matrix_identity.f90)
target_link_libraries(matrix_identity
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_matrix)

# Test with the uniform matrix
set_source_files_properties(matrix_uniform.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_uniform
    matrix_uniform.f90)
target_link_libraries(matrix_uniform
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_matrix)

# Test with the banded uniform matrix
set_source_files_properties(matrix_banded_uniform.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(matrix_banded_uniform
    matrix_banded_uniform.f90)
target_link_libraries(matrix_banded_uniform
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_matrix)

# Test with random rounded tensor
set_source_files_properties(tensor_integer.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(tensor_integer
    tensor_integer.f90)
target_link_libraries(tensor_integer
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_tensor)

# Test with eye tensor
set_source_files_properties(tensor_eye.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(tensor_eye
    tensor_eye.f90)
target_link_libraries(tensor_eye
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_tensor)

# Test with uniform tensor
set_source_files_properties(tensor_uniform.f90
    PROPERTIES Fortran_PREPROCESS ON)
add_executable(tensor_uniform
    tensor_uniform.f90)
target_link_libraries(tensor_uniform
    PRIVATE maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_comparison
    maria_la_utils maria_la_core maximum_tensor)
