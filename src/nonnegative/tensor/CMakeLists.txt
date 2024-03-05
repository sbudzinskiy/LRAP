add_library(nonnegative_tensor
    STATIC nonnegative_tensor_mod.f90 ./sub/nonnegative_tensor_sub.f90)
target_link_libraries(nonnegative_tensor
    PRIVATE maria_kinds maria_constants maria_comparison maria_reports 
    maria_la_core maria_lr_tsvd maria_lr_cross maria_lr_maxvol maria_tt_utils maria_tt_tsvd)
