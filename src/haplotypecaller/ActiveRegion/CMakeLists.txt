set(PACK_NAME region)
set(PACK_SRC ActiveMainThread.cpp ActiveBase.cpp HcActiveBase.cpp reference_manager.cpp ActiveRegionEngine.cpp RegionManager.cpp dbsnp_manager.cpp)
set(PACK_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR} ../common ../common/assemble ../common/utils/ ../common/base/include ../../interface
                     ../../common ../../IOStream ${Boost_INCLUDE_DIRS} ${SPDLOG_INCLUDE_DIR} ../../rovaca_logger ../bqsr)
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB hts utils IOStream logger)

add_library(${PACK_NAME} SHARED ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

if(ENABLE_GOOGLE_TEST)
    set(test_name ActiveRegionTest)
    set(test_src ActiveRegionTest.cpp)
    add_executable(${test_name} ${test_src})
    target_include_directories(${test_name} PUBLIC ${PACK_INCLUDE_DIR})
    target_link_libraries(${test_name} hts utils IOStream region)

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
