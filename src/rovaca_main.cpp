#include "rovaca_main.h"

#include <iostream>

#include "rovaca_logger.h"
#include "rovaca_signal_handler.h"
#include "rovaca_tool.hpp"

int rovaca_main(int argc, char* argv[])
{
    RovacaLogger::init_assemble_ptr();
    RovacaLogger::set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] [%s:%!:%#] %v.");

    auto& register_ = RovacaToolRegister::instance();
    if (argc < 2) {
        RovacaLogger::error("no rovaca tool specified");
        return register_.supported_tools();
    }

    std::unique_ptr<RovacaTool> run_tool = register_.creat_tool(argv[1]);
    if (!run_tool) {
        RovacaLogger::error("invalid rovaca tool name: {}", argv[1]);
        return register_.supported_tools();
    }

    if (!run_tool->initialize_args(argc, argv)) {
        RovacaLogger::error("invalid rovaca tool args");
        return EXIT_FAILURE;
    }

    try {
        RovacaSignalHandler::setup_signal_handler();
        run_tool->run();
    }
    catch (const RovacaSignalException& e) {
        RovacaLogger::error("{}", e.what());
        run_tool->clear_and_exit();
    }
    catch (const std::runtime_error& e) {
        RovacaLogger::error("{}", e.what());
        run_tool->clear_and_exit();
    }

    return EXIT_SUCCESS;
}