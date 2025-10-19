#ifndef Haplotypecaller_H
#define Haplotypecaller_H

#include <memory>
#include <thread>

#include "BlockingQueue.h"
#include "bqsr/bqsr_read_transformer.h"
#include "haplotypecaller_args.h"
#include "haplotypecaller_engine.h"
#include "rovaca_tool.hpp"
#include "writer/writer.h"
class HaplotypeCaller : public RovacaTool
{
private:
    std::unique_ptr<HaplotypeCallerArgs> hc_args;
    std::unique_ptr<HaplotypeCallerEngine> hc_engine;
    std::unique_ptr<std::thread> assemble_output_thread_;
    std::unique_ptr<rovaca::Writer> writer_;
    BQSRReadTransformer* apply_bqsr_;
    bool debug_assemble_stream_;

public:
    HaplotypeCaller();
    ~HaplotypeCaller();
    void do_work() override;
    void clear_and_exit() override;

private:
    const UniqueStream& custom_streamer() override;

    void init_args();
};

#endif  // Haplotypecaller_H