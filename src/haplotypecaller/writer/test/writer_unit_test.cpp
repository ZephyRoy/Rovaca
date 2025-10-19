#include <gtest/gtest.h>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

#include "common/utils/BlockingQueue.h"
#include "genotype/forward.h"
#include "genotype/genotype_argument.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "writer.h"

using namespace rovaca;
using namespace std;

const char* test_output = "/data/rovaca-dev/zhangtiefeng/dev/hc/build/writer_test.g.vcf.gz";
const char* fasta_file = "/data/rovaca-dev/sungeng/Datas/hg19.fasta";
const char* bam_file = "/data/rovaca-dev/sungeng/Datas/NA12878.sort.dup.bam";

static std::vector<pWriterTask> create_task(Writer* writer, const shared_ptr<char>& ref)
{
    /*
    01 del 跨越下一个source
    02 del 跨越两个source
    04 del 跨越两个source且pos位数增加(三位变四位)：999->1050
    */

    // 期待: 1
    const char* s1 = "chrM\t1\t.\tCCCCC\tT,<NON_REF>\t275\t.\tDP=6658\tGT\t1/1\n";

    // 期待: 6-72
    const char* s2 = "chrM\t3\t.\tG\t<NON_REF>\t.\t.\tEND=72\tGT\t0/0\n";

    // 期待: 73, 83-92
    const char* s3 = "chrM\t73\t.\tATGCATGCTG\tT,<NON_REF>\t275\t.\tDP=6658\tGT\t1/1\n";
    const char* s4 = "chrM\t75\t.\tG\t<NON_REF>\t.\t.\tEND=77\tGT\t0/0\n";
    const char* s5 = "chrM\t78\t.\tG\t<NON_REF>\t.\t.\tEND=80\tGT\t0/0\n";
    const char* s6 = "chrM\t81\t.\tG\t<NON_REF>\t.\t.\tEND=92\tGT\t0/0\n";

    // 期待: 93, 103-120
    const char* s7 = "chrM\t93\t.\tATGCATGCTG\tT,<NON_REF>\t275\t.\tDP=6658\tGT\t1/1\n";
    const char* s8 = "chrM\t94\t.\tG\t<NON_REF>\t.\t.\tEND=96\tGT\t0/0\n";
    const char* s9 = "chrM\t97\t.\tG\t<NON_REF>\t.\t.\tEND=98\tGT\t0/0\n";
    const char* s10 = "chrM\t99\t.\tG\t<NON_REF>\t.\t.\tEND=120\tGT\t0/0\n";

    std::vector<pWriterTask> result;

    int32_t tid = 0, source_id = 0;

    kstring_t* cache;
    pWriterTask task;

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s1);
    memcpy(cache->s, s1, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 1, 6, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s2);
    memcpy(cache->s, s2, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s3);
    memcpy(cache->s, s3, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 1, 83, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s4);
    memcpy(cache->s, s4, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s5);
    memcpy(cache->s, s5, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s6);
    memcpy(cache->s, s6, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s7);
    memcpy(cache->s, s7, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 1, 103, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s8);
    memcpy(cache->s, s8, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s9);
    memcpy(cache->s, s9, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    cache = writer->pop_bcf_cache();
    cache->l = strlen(s10);
    memcpy(cache->s, s10, cache->l * sizeof(char));
    task = new WriterTask{tid, source_id++, cache, 0, -1, ref};
    result.push_back(task);

    return result;
}

TEST(WriterUnitTest, resolveOverridesAndWriteTask)
{
    samFile* bam = sam_open(bam_file, "r");
    if (nullptr == bam) {
        printf("sam_open: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    bam_hdr_t* hdr = sam_hdr_read(bam);
    if (nullptr == hdr) {
        printf("sam_hdr_read: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    faidx_t* fai = fai_load(fasta_file);
    if (nullptr == fai) {
        printf("fai_load: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    int32_t len;
    shared_ptr<char> ref(faidx_fetch_seq(fai, "chrM", 0, 16000, &len));
    if (nullptr == ref) {
        printf("faidx_fetch_seq: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    pHCArgs args = new GenotypeArgument{};
    args->output = test_output;
    args->tool_name = "HaplotypeCaller";
    args->command_line = "rovaca HaplotypeCaller WriterUnitTest";
    for (int32_t i = 1; i <= 60; ++i) {
        args->gvcf_gq_bands.push_back(i);
    }
    args->gvcf_gq_bands.push_back(70);
    args->gvcf_gq_bands.push_back(80);
    args->gvcf_gq_bands.push_back(90);
    args->gvcf_gq_bands.push_back(90);
    args->init_reference_confidence_mode(ReferenceConfidenceMode::GVCF);

    htsThreadPool pool{nullptr, 0};
    pool.pool = hts_tpool_init(4);

    BlockingQueue<pWriterTask>* result_queue = new BlockingQueue<pWriterTask>{256};

    Writer writer{false, 6, args, hdr, &pool, result_queue};
    std::vector<pWriterTask> tasks = create_task(&writer, ref);
    for_each(tasks.begin(), tasks.end(), [=](pWriterTask t) { result_queue->push(t); });

    writer.start();
    writer.update_last_id_back(int32_t(tasks.size()));
    writer.stop();

    delete result_queue;
    hts_tpool_destroy(pool.pool);
    delete args;
    fai_destroy(fai);
    sam_hdr_destroy(hdr);
    sam_close(bam);
}