add_library(nn_matrix_completion
    STATIC nn_matrix_completion_mod.f90 ./sub/nn_matrix_completion_sub.f90)
target_link_libraries(nn_matrix_completion
    PRIVATE maria_kinds maria_constants maria_comparison maria_reports 
    maria_la_core maria_lr_tsvd maria_lr_cross maria_lr_maxvol maria_lr_geom)
