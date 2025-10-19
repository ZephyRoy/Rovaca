#ifndef ROVACA_HC_ROVACA_LOGGER_H_
#define ROVACA_HC_ROVACA_LOGGER_H_
#include "spdlog/spdlog.h"

/*!
 * @brief Adapt the log function for the assemble part, pass the function pointer to assemble. Variadic arguments are not supported here, only for printing simple messages
 * @note Need to use macros to get file name, function name, and line number. The macro is defined in src/haplotypecaller/common/base/include/debug.h
 */
extern "C"
{
    void rovacalog_trace(int line, const char* file_name, const char* func_name, const char* msg);
    void rovacalog_debug(int line, const char* file_name, const char* func_name, const char* msg);
    void rovacalog_info(int line, const char* file_name, const char* func_name, const char* msg);
    void rovacalog_warn(int line, const char* file_name, const char* func_name, const char* msg);
    void rovacalog_error(int line, const char* file_name, const char* func_name, const char* msg);
    void rovacalog_critical(int line, const char* file_name, const char* func_name, const char* msg);
}

struct RovacaLocation
{
    /*!
     * @brief Use GCC internal instructions to get file name, function name, and line number
     * @param line
     * @param file_name
     * @param func_name
     */
    constexpr explicit RovacaLocation(int line = __builtin_LINE(), const char* file_name = __builtin_FILE(),
                                    const char* func_name = __builtin_FUNCTION()) noexcept
        : line_(line)
        , file_name_(file_name)
        , function_name_(func_name)
    {}

    /*!
     * @brief RovacaLocation changed spdlog::source_loc，easy to use spdlog print function related information
     * @return
     */
    [[nodiscard]] spdlog::source_loc to_spdlog_source() const
    {
        const char* file_name_no_path = get_file_name_no_path(file_name_);
        return {file_name_no_path, line_, function_name_};
    }

    /*!
     * @brief Get the file name without path
     * @param file_name
     * @return
     */
    static const char* get_file_name_no_path(const char* file_name)
    {
        const char* file_name_no_path = strrchr(file_name, '/');
        file_name_no_path = file_name_no_path == nullptr ? file_name : file_name_no_path + 1;
        return file_name_no_path;
    }

    int line_;
    const char* file_name_;
    const char* function_name_;
};

/*!
 * @brief Implemented using C++17 deduction guide method: variadic + default parameters
 */
namespace RovacaLogger
{

void init_assemble_ptr();
void set_level(spdlog::level::level_enum log_level);
void set_pattern(std::string pattern, spdlog::pattern_time_type time_type = spdlog::pattern_time_type::local);

template <typename... Args>
struct trace
{
    constexpr explicit trace(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::trace, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
trace(fmt::format_string<Args...> fmt, Args&&... args) -> trace<Args...>;

template <typename... Args>
struct debug
{
    constexpr explicit debug(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::debug, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
debug(fmt::format_string<Args...> fmt, Args&&... args) -> debug<Args...>;

template <typename... Args>
struct info
{
    constexpr explicit info(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::info, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
info(fmt::format_string<Args...> fmt, Args&&... args) -> info<Args...>;

template <typename... Args>
struct warn
{
    constexpr explicit warn(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::warn, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
warn(fmt::format_string<Args...> fmt, Args&&... args) -> warn<Args...>;

template <typename... Args>
struct error
{
    constexpr explicit error(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::err, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
error(fmt::format_string<Args...> fmt, Args&&... args) -> error<Args...>;

template <typename... Args>
struct critical
{
    constexpr explicit critical(fmt::format_string<Args...> fmt, Args&&... args, RovacaLocation location = RovacaLocation{})
    {
        spdlog::log(location.to_spdlog_source(), spdlog::level::critical, fmt, std::forward<Args>(args)...);
    }
};

template <typename... Args>
critical(fmt::format_string<Args...> fmt, Args&&... args) -> critical<Args...>;

}  // namespace RovacaLogger

/**
 * @brief Check if the condition is met, if so, the program exits - stricter
 * @param cond bool type condition
 * @param mesg const char* error message
 */
#define CHECK_CONDITION_EXIT(cond, fmt, ...)       \
    do {                                           \
        if (__glibc_unlikely(cond)) {              \
            RovacaLogger::error(fmt, ##__VA_ARGS__); \
            exit(EXIT_FAILURE);                    \
        }                                          \
    } while (false)

#endif  // ROVACA_HC_ROVACA_LOGGER_H_
