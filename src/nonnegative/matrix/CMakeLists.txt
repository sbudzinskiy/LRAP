add_library(nonnegative_matrix
    STATIC nonnegative_matrix_mod.f90 ./sub/nonnegative_matrix_sub.f90)
target_link_libraries(nonnegative_matrix
    PRIVATE maria_kinds maria_constants maria_comparison maria_reports 
    maria_la_core maria_lr_tsvd maria_lr_cross maria_lr_maxvol)
