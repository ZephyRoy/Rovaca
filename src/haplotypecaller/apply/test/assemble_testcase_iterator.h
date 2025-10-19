#ifndef ASSEMBLE_REGION_TESTCASE_ITERATOR_H
#define ASSEMBLE_REGION_TESTCASE_ITERATOR_H

#include <cmath>
#include <string>
#include <vector>

#include "htslib/sam.h"
#include "testcase_iterator.h"

static const std::string hdr =
    "data:,"
    "@CO\tThis is for testing.\n"
    "@SQ\tSN:chrM\tLN:16571\n"
    "@SQ\tSN:chr1\tLN:249250621\n"
    "@SQ\tSN:chr2\tLN:243199373\n"
    "@SQ\tSN:chr3\tLN:198022430\n"
    "@SQ\tSN:chr4\tLN:191154276\n"
    "@SQ\tSN:chr5\tLN:180915260\n"
    "@SQ\tSN:chr6\tLN:171115067\n"
    "@SQ\tSN:chr7\tLN:159138663\n"
    "@SQ\tSN:chr8\tLN:146364022\n"
    "@SQ\tSN:chr9\tLN:141213431\n"
    "@SQ\tSN:chr10\tLN:135534747\n"
    "@SQ\tSN:chr11\tLN:135006516\n"
    "@SQ\tSN:chr12\tLN:133851895\n"
    "@SQ\tSN:chr13\tLN:115169878\n"
    "@SQ\tSN:chr14\tLN:107349540\n"
    "@SQ\tSN:chr15\tLN:102531392\n"
    "@SQ\tSN:chr16\tLN:90354753\n"
    "@SQ\tSN:chr17\tLN:81195210\n"
    "@SQ\tSN:chr18\tLN:78077248\n"
    "@SQ\tSN:chr19\tLN:59128983\n"
    "@SQ\tSN:chr20\tLN:63025520\n"
    "@SQ\tSN:chr21\tLN:48129895\n"
    "@SQ\tSN:chr22\tLN:51304566\n"
    "@SQ\tSN:chrX\tLN:155270560\n"
    "@SQ\tSN:chrY\tLN:59373566\n"
    "@SQ\tSN:chr1_gl000191_random\tLN:106433\n"
    "@SQ\tSN:chr1_gl000192_random\tLN:547496\n"
    "@SQ\tSN:chr4_ctg9_hap1\tLN:590426\n"
    "@SQ\tSN:chr4_gl000193_random\tLN:189789\n"
    "@SQ\tSN:chr4_gl000194_random\tLN:191469\n"
    "@SQ\tSN:chr6_apd_hap1\tLN:4622290\n"
    "@SQ\tSN:chr6_cox_hap2\tLN:4795371\n"
    "@SQ\tSN:chr6_dbb_hap3\tLN:4610396\n"
    "@SQ\tSN:chr6_mann_hap4\tLN:4683263\n"
    "@SQ\tSN:chr6_mcf_hap5\tLN:4833398\n"
    "@SQ\tSN:chr6_qbl_hap6\tLN:4611984\n"
    "@SQ\tSN:chr6_ssto_hap7\tLN:4928567\n"
    "@SQ\tSN:chr7_gl000195_random\tLN:182896\n"
    "@SQ\tSN:chr8_gl000196_random\tLN:38914\n"
    "@SQ\tSN:chr8_gl000197_random\tLN:37175\n"
    "@SQ\tSN:chr9_gl000198_random\tLN:90085\n"
    "@SQ\tSN:chr9_gl000199_random\tLN:169874\n"
    "@SQ\tSN:chr9_gl000200_random\tLN:187035\n"
    "@SQ\tSN:chr9_gl000201_random\tLN:36148\n"
    "@SQ\tSN:chr11_gl000202_random\tLN:40103\n"
    "@SQ\tSN:chr17_ctg5_hap1\tLN:1680828\n"
    "@SQ\tSN:chr17_gl000203_random\tLN:37498\n"
    "@SQ\tSN:chr17_gl000204_random\tLN:81310\n"
    "@SQ\tSN:chr17_gl000205_random\tLN:174588\n"
    "@SQ\tSN:chr17_gl000206_random\tLN:41001\n"
    "@SQ\tSN:chr18_gl000207_random\tLN:4262\n"
    "@SQ\tSN:chr19_gl000208_random\tLN:92689\n"
    "@SQ\tSN:chr19_gl000209_random\tLN:159169\n"
    "@SQ\tSN:chr21_gl000210_random\tLN:27682\n"
    "@SQ\tSN:chrUn_gl000211\tLN:166566\n"
    "@SQ\tSN:chrUn_gl000212\tLN:186858\n"
    "@SQ\tSN:chrUn_gl000213\tLN:164239\n"
    "@SQ\tSN:chrUn_gl000214\tLN:137718\n"
    "@SQ\tSN:chrUn_gl000215\tLN:172545\n"
    "@SQ\tSN:chrUn_gl000216\tLN:172294\n"
    "@SQ\tSN:chrUn_gl000217\tLN:172149\n"
    "@SQ\tSN:chrUn_gl000218\tLN:161147\n"
    "@SQ\tSN:chrUn_gl000219\tLN:179198\n"
    "@SQ\tSN:chrUn_gl000220\tLN:161802\n"
    "@SQ\tSN:chrUn_gl000221\tLN:155397\n"
    "@SQ\tSN:chrUn_gl000222\tLN:186861\n"
    "@SQ\tSN:chrUn_gl000223\tLN:180455\n"
    "@SQ\tSN:chrUn_gl000224\tLN:179693\n"
    "@SQ\tSN:chrUn_gl000225\tLN:211173\n"
    "@SQ\tSN:chrUn_gl000226\tLN:15008\n"
    "@SQ\tSN:chrUn_gl000227\tLN:128374\n"
    "@SQ\tSN:chrUn_gl000228\tLN:129120\n"
    "@SQ\tSN:chrUn_gl000229\tLN:19913\n"
    "@SQ\tSN:chrUn_gl000230\tLN:43691\n"
    "@SQ\tSN:chrUn_gl000231\tLN:27386\n"
    "@SQ\tSN:chrUn_gl000232\tLN:40652\n"
    "@SQ\tSN:chrUn_gl000233\tLN:45941\n"
    "@SQ\tSN:chrUn_gl000234\tLN:40531\n"
    "@SQ\tSN:chrUn_gl000235\tLN:34474\n"
    "@SQ\tSN:chrUn_gl000236\tLN:41934\n"
    "@SQ\tSN:chrUn_gl000237\tLN:45867\n"
    "@SQ\tSN:chrUn_gl000238\tLN:39939\n"
    "@SQ\tSN:chrUn_gl000239\tLN:33824\n"
    "@SQ\tSN:chrUn_gl000240\tLN:41933\n"
    "@SQ\tSN:chrUn_gl000241\tLN:42152\n"
    "@SQ\tSN:chrUn_gl000242\tLN:43523\n"
    "@SQ\tSN:chrUn_gl000243\tLN:43341\n"
    "@SQ\tSN:chrUn_gl000244\tLN:39929\n"
    "@SQ\tSN:chrUn_gl000245\tLN:36651\n"
    "@SQ\tSN:chrUn_gl000246\tLN:38154\n"
    "@SQ\tSN:chrUn_gl000247\tLN:36422\n"
    "@SQ\tSN:chrUn_gl000248\tLN:39786\n"
    "@SQ\tSN:chrUn_gl000249\tLN:38502\n"
    "@RG\tID:NA12878\tPL:COMPLETE\tPU:NA12878\tLB:NA12878\tSM:SF3\tCN:BGI\n";

class AssembleRegionCase
{
public:
    std::string contig;
    hts_pos_t region_beg;
    hts_pos_t region_end;
    int region_size;
    std::vector<bam1_t*> region_reads;
    std::string reference;
    AssembleRegionCase()
        : contig("")
        , region_beg(0)
        , region_end(0)
        , region_size(0)
        , region_reads({})
        , reference("")
    {}
    AssembleRegionCase(const std::string& chr, hts_pos_t beg, hts_pos_t end, int size, const std::vector<bam1_t*>& reads,
                       const std::string& ref)
        : contig(chr)
        , region_beg(beg)
        , region_end(end)
        , region_size(size)
        , region_reads(reads)
        , reference(ref)
    {}
};

class AssembleTestCaseIterator : public TestCaseIterator<AssembleRegionCase>
{
private:
    samFile* fp;
    bam_hdr_t* hp;

public:
    AssembleTestCaseIterator(std::istream* const input)
        : TestCaseIterator<AssembleRegionCase>(input)
        , fp(sam_open(hdr.c_str(), "r"))
        , hp(sam_hdr_read(fp))
    {
        current = fetch_next();
    }
    ~AssembleTestCaseIterator()
    {
        // 每次构造fp->line时，.s内存为user define,因此默认析构前malloc,防止double free.
        fp->line.s = (char*)malloc(1);
        bam_hdr_destroy(hp);
        sam_close(fp);
    }
    AssembleRegionCase fetch_next() override;
};

AssembleRegionCase AssembleTestCaseIterator::fetch_next()
{
    if (input_stream == nullptr) {
        return AssembleRegionCase{};
    }

    std::string read_line;
    std::getline(*input_stream, read_line);
    if (read_line.empty()) {
        input_stream = nullptr;
        std::cerr << "All " << case_count << " cases have been loaded." << std::endl;
        return AssembleRegionCase{};
    }
    accumulate_cases();
    int64_t beg, end;
    int size;
    char contig[32];
    sscanf(read_line.c_str(), "%s\t%ld\t%ld\t%d", contig, &beg, &end, &size);
    std::vector<bam1_t*> reads{};

    // 这里使用sam_read1()::hts_get_line函数,在解码类似百分号（%）后跟两个十六进制数字这种字符串时出现错误。
    // 例如在URL编码中，空格字符编码为 %20，而字符 u 的编码为 %75。
    int count = size;
    while (count--) {
        std::getline(*input_stream, read_line);
        bam1_t* bam = bam_init1();
        fp->line = {0, 0, const_cast<char*>(read_line.c_str())};
        if (sam_parse1(&fp->line, hp, bam) < 0) {
            std::cerr << "get read error!" << std::endl;
        }
        reads.emplace_back(bam);
    }
    std::getline(*input_stream, read_line);
    AssembleRegionCase testcase{contig, beg, end, size, reads, read_line};
    return testcase;
}
#endif  // ASSEMBLE_REGION_TESTCASE_ITERATOR_H