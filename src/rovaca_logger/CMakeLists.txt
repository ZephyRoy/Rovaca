set(PACK_NAME logger)
set(PACK_SRC rovaca_logger.cpp rovaca_logger.h)
set(PACK_INCLUDE_DIR ${SPDLOG_INCLUDE_DIR})
set(PACK_LINK_DIR ${SPDLOG_LIBRARY_DIR})
set(PACK_LINK_LIB spdlog)

add_library(${PACK_NAME} SHARED ${PACK_SRC})
target_include_directories(${PACK_NAME} PUBLIC ${PACK_INCLUDE_DIR})
target_link_directories(${PACK_NAME} PUBLIC ${PACK_LINK_DIR})
target_link_libraries(${PACK_NAME} PUBLIC ${PACK_LINK_LIB})

install(TARGETS ${PACK_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
