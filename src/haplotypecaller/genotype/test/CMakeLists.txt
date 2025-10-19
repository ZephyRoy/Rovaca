file(GLOB PACK_FILES *.cpp)

set(PACK_NAME genotype_test)
set(PACK_SRC ${PACK_FILES})
set(PACK_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR} .. ../../common/utils ../../common/base/include ../../common/assemble ../../../rovaca_logger ${SPDLOG_INCLUDE_DIR})
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB hts gtest gtest_main genotype utils writer base logger)

foreach(TEST_FILE ${PACK_FILES})
    get_filename_component(FILE_NAME ${TEST_FILE} NAME_WE)

    set(TEST_NAME ${PACK_NAME}_${FILE_NAME})

    add_executable(${TEST_NAME} ${TEST_FILE})
    target_include_directories(${TEST_NAME} PRIVATE ${PACK_INCLUDE_DIR})
    target_link_libraries(${TEST_NAME} PRIVATE ${PACK_LINK_LIB})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
