set(PACK_NAME active_region_test)
set(PACK_SRC TestActiveRegionEngine.cpp TestHCActiveBase.cpp TestActiveMainThread.cpp TestReferenceManager.cc)
set(PACK_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR} .. ../../IOStream)
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB hts gtest_main gtest utils region IOStream)

add_executable(${PACK_NAME} ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

add_test(NAME ${PACK_NAME} COMMAND ${PACK_NAME})
