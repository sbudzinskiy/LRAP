add_library(maximum_matrix
    STATIC maximum_matrix_mod.f90 ./sub/maximum_matrix_sub.f90)
target_link_libraries(maximum_matrix
    PRIVATE maria_kinds maria_constants maria_comparison maria_reports 
    maria_la_core maria_la_utils maria_lr_tsvd maria_lr_cross maria_lr_maxvol)
