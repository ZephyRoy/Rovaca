set(PACK_NAME utils)
set(PACK_SRC math_utils.cpp math_utils.h quality_utils.cpp quality_utils.h bam_data_pool.hpp rovaca_memory_pool.cpp rovaca_memory_pool.h)
set(PACK_INCLUDE_DIR ../base/include ${HTSLIB_INCLUDE_DIR} ${SPDLOG_INCLUDE_DIR} ../../../rovaca_logger)
set(PACK_LINK_DIR)
set(PACK_LINK_LIB pthread logger)

add_library(${PACK_NAME} SHARED ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

if(ENABLE_GOOGLE_TEST)
    add_subdirectory(test)
endif()

if(ENABLE_ASAN)
    include(asan)
    enable_sanitizer(${PACK_NAME})
endif()

if(ENABLE_GCOV)
    include(gcov)
    enable_gcov(${PACK_NAME})
endif()

install(TARGETS ${PACK_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin)
