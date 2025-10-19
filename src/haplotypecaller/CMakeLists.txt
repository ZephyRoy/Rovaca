add_subdirectory(genotype)
add_subdirectory(common)
add_subdirectory(apply)
add_subdirectory(ActiveRegion)
add_subdirectory(pairhmm)
add_subdirectory(writer)
add_subdirectory(bqsr)

set(PACK_NAME haplotypecaller)
set(PACK_SRC haplotypecaller.cpp haplotypecaller_engine.cpp haplotypecaller_args.cpp)
set(PACK_INCLUDE_DIR .. common ../interface ../IOStream ../common common/base/include common/assemble common/utils/ ActiveRegion genotype
                     apply pairhmm apply/common ${Boost_INCLUDE_DIRS} ${HTSLIB_INCLUDE_DIR})
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB genotype pairhmm region utils base assemble writer bqsr)

add_library(${PACK_NAME} SHARED ${PACK_SRC})
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
