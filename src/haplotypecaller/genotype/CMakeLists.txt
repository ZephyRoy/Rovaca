file(GLOB UTILS_SRC "utils/*.cpp")
file(GLOB INTERFACE_SRC "interface/*.hpp")
file(GLOB ANNOTATION_INFO_SRC "annotation/info/*.cpp")
file(GLOB ANNOTATION_FORMAT_SRC "annotation/format/*.cpp")
file(GLOB ALL_SRC "*.cpp")

set(PACK_NAME genotype)
set(PACK_SRC ${UTILS_SRC} ${INTERFACE_SRC} ${ANNOTATION_INFO_SRC} ${ANNOTATION_FORMAT_SRC} ${ALL_SRC})
set(PACK_INCLUDE_DIR ${HTSLIB_INCLUDE_DIR} ../ ../../ ../common/utils ../common/assemble ${SPDLOG_INCLUDE_DIR} ../../rovaca_logger)
set(PACK_LINK_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(PACK_LINK_LIB hts utils base logger)

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
