set(PACK_NAME IOStream_test)
set(PACK_SRC bam_loader_test.cc bed_loader_test.cc fasta_loader_test.cc ../bam_loader.cpp ../fasta_loader.cpp)
set(PACK_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR} ${SPDLOG_INCLUDE_DIR} ..)
set(PACK_LINK_DIR ${HTSLIB_LIBRARY_DIR})
set(PACK_LINK_LIB hts gtest_main gtest pthread)

add_executable(${PACK_NAME} ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

add_test(NAME ${PACK_NAME} COMMAND ${PACK_NAME})
