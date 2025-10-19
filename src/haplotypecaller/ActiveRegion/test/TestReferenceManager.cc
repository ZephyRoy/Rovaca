#include <string>

#include "gtest/gtest.h"
#include "reference_manager.h"

TEST(ReferenceManager, case1)
{
    std::string ref("/data/pipelines/WGS_bgionline/db/GRCh37/ref/GRCh37_no_alt.fna");
    ReferenceManager inst(ref);
    std::thread thread_manager(&ReferenceManager::run, &inst);

    int length = 0;
    for (int i = 0; i < 85; i++) {
        std::shared_ptr<char> base = inst.get(i, length);
        inst.pop();
        EXPECT_NE(base, nullptr);
    }
    thread_manager.join();
}

TEST(ReferenceManager, case2)
{
    std::string ref("/data/pipelines/WGS_bgionline/db/GRCh37/ref/GRCh37_no_alt.fna");
    ReferenceManager inst(ref);
    std::thread thread_manager(&ReferenceManager::run, &inst);

    int length = 0;
    for (int i = 0; i < 2; i++) {
        std::shared_ptr<char> base = inst.get(i, length);
        EXPECT_NE(base, nullptr);
    }
    inst.assign_finish();
    thread_manager.join();
}