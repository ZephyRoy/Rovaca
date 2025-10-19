set(PACK_NAME writer_unit_test)
set(PACK_SRC writer_unit_test.cpp)
set(PACK_INCLUDE_DIR .. ../.. ${HTSLIB_INCLUDE_DIR})
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB gtest_main gtest pthread writer genotype)

add_executable(${PACK_NAME} ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

add_test(NAME ${PACK_NAME} COMMAND ${PACK_NAME})
