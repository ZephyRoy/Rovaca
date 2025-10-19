set(PACK_NAME base)
set(PACK_SRC lib/mempool/mem_pool_fast.c lib/mempool/mem_pool_auto_enlarge.c lib/mempool/ring_mem_pool.c lib/sw_avx/sw_avx2.c
             lib/sw_avx/sw_avx.c)
set(PACK_INCLUDE_DIR include ${HTSLIB_INCLUDE_DIR})
set(PACK_LINK_DIR)
set(PACK_LINK_LIB pthread m)

add_library(${PACK_NAME} SHARED ${PACK_SRC})
target_compile_options(${PACK_NAME} PRIVATE -mavx -mavx2 -fgnu89-inline)
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

if(ENABLE_ASAN)
    include(asan)
    enable_sanitizer(${PACK_NAME})
endif()

if(ENABLE_GCOV)
    include(gcov)
    enable_gcov(${PACK_NAME})
endif()

install(TARGETS ${PACK_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin)
