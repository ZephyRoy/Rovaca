#ifndef ROVACA_SIGNAL_EXCEPTHON_H
#define ROVACA_SIGNAL_EXCEPTHON_H

#include <csignal>
#include <stdexcept>

#include "rovaca_logger.h"
#include "rovaca_tool.hpp"

class RovacaSignalException : public std::runtime_error
{
public:
    explicit RovacaSignalException(const std::string& message)
        : std::runtime_error(message)
    {}
};

class RovacaSignalHandler
{
private:
    RovacaSignalHandler() = default;

public:
    static void setup_signal_handler();
    static void signal_handler(int signal);
};

void RovacaSignalHandler::setup_signal_handler()
{
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
}

void RovacaSignalHandler::signal_handler(int signal)
{
    if (signal == SIGINT) {
        RovacaLogger::error("recieved signal: SIGINT. exiting");
        if (k_outfile != nullptr) hts_close(k_outfile);
        _exit(-1);
    }
    if (signal == SIGTERM) {
        RovacaLogger::error("recieved signal: SIGTERM. exiting");
        if (k_outfile != nullptr) hts_close(k_outfile);
        _exit(-1);
    }
}

#endif  // Rovaca_SIGNAL_EXCEPTHON_H