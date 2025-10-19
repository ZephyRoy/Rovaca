set(PackName bqsr)
set(PackSrc
    bqsr_read_covarivates.cpp
    bqsr_read_transformer.cpp
    apply_bqsr_utils.cpp
)

add_library(
    ${PackName}
    SHARED
    ${PackSrc}
)

target_link_libraries(
    ${PackName}
    PUBLIC
)

target_include_directories(
    ${PackName}
    PUBLIC
    ..
    ../common/
    ../../IOStream/
    ../../rovaca_logger
    ../common/utils/
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIR}
    ${SPDLOG_INCLUDE_DIR})
    
    install(TARGETS ${PackName} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin)