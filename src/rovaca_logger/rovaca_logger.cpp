#include "rovaca_logger.h"

/*!
 * @brief The function pointer here is used for assemble. Export these variables in the common assemble header file, then link this rovaca_logger.so
 */
typedef void (*p_rovaca_log_func)(int, const char *, const char *, const char *);
p_rovaca_log_func rovaca_assemble_log_trace;
p_rovaca_log_func rovaca_assemble_log_debug;
p_rovaca_log_func rovaca_assemble_log_info;
p_rovaca_log_func rovaca_assemble_log_warn;
p_rovaca_log_func rovaca_assemble_log_error;
p_rovaca_log_func rovaca_assemble_log_critical;

void rovacalog_trace(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::trace, msg);
}

void rovacalog_debug(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::debug, msg);
}

void rovacalog_info(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::info, msg);
}

void rovacalog_warn(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::warn, msg);
}

void rovacalog_error(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::err, msg);
}

void rovacalog_critical(int line, const char *file_name, const char *func_name, const char *msg)
{
    spdlog::source_loc loc = RovacaLocation{line, file_name, func_name}.to_spdlog_source();
    spdlog::log(loc, spdlog::level::critical, msg);
}

void RovacaLogger::init_assemble_ptr()
{
    rovaca_assemble_log_trace = &rovacalog_trace;
    rovaca_assemble_log_debug = &rovacalog_debug;
    rovaca_assemble_log_info = &rovacalog_info;
    rovaca_assemble_log_warn = &rovacalog_warn;
    rovaca_assemble_log_error = &rovacalog_error;
    rovaca_assemble_log_critical = &rovacalog_critical;
}

void RovacaLogger::set_level(spdlog::level::level_enum log_level) { spdlog::set_level(log_level); }

void RovacaLogger::set_pattern(std::string pattern, spdlog::pattern_time_type time_type)
{
    spdlog::set_pattern(std::move(pattern), time_type);
}
