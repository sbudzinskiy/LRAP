add_library(maximum_tensor
    STATIC maximum_tensor_mod.f90 ./sub/maximum_tensor_sub.f90)
target_link_libraries(maximum_tensor
    PRIVATE maria_kinds maria_constants maria_comparison maria_reports 
    maria_la_core maria_tt_tsvd maria_tt_utils)
