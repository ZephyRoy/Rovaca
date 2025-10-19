#ifndef RUN_ASSEMBLER_H
#define RUN_ASSEMBLER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "assemble_argument.h"
#include "assemble_engine.h"
#include "assemble_testcase_iterator.h"
#include "hc_assemble_main.h"
#include "htslib/sam.h"
#include "testcase_iterator.h"
#include "testcase_loader.h"

class FakeAssembler
{
private:
    AssembleArgument arguments_;
    std::fstream assemble_debug_streamer_;
    TestCaseLoader<AssembleTestCaseIterator> test_case_;
    pMemoryPool target_mem_;
    std::string chr_ref;
    uint32_t chr_ref_len;

public:
    FakeAssembler() = delete;
    FakeAssembler(const char* casePath, const char* refPath, const char* resultPath, const AssembleArgument& argument,
                  pMemoryPool targetMem);
    ~FakeAssembler();
    void run();
};

#endif  // RUN_ASSEMBLER_H