#ifndef ROVACA_HC_DEBUG_UTILS_H_
#define ROVACA_HC_DEBUG_UTILS_H_
#include "../forward.h"

namespace rovaca
{

/*!
 * @brief 集成测试时，对象 <-> 字符串 互相转换
 */
namespace DebugUtils
{

std::string cigar2str(pCigar c);
std::string cigar2str(const uint32_t* cigar, uint32_t num);
pCigar str2cigar(const std::string& s, pMemoryPool pool);

std::string bases2str(pBases b);
std::string bases2str(const uint8_t* b, int32_t len);
pBases str2bases(const std::string& s, pMemoryPool pool);

/*!
 * @brief str 格式：chr1:1-100
 */
std::string interval2str(pInterfaceLocatable i, bam_hdr_t* header);
pSimpleInterval str2interval(const std::string& s, bam_hdr_t* header, pMemoryPool pool);

/*!
 * @brief 单倍体格式：h cigar interval score ali kmer is_ref bases
 * score可能为nan
 */
std::string haplotype2str(pHaplotype h, bam_hdr_t* header);
pHaplotype str2haplotype(const std::string& s, bam_hdr_t* header, pMemoryPool pool);

/*!
 * @brief Read格式：name flag interval mq cigar mtid mpos isize seq qual
 * 非成对read，
 */
std::string read2str(pReadRecord r, bam_hdr_t* header, pMemoryPool pool);
pReadRecord str2read(const std::string& s, bam_hdr_t* header, pMemoryPool pool, pBamDataPool bpool);

std::string variant2str(pVariant v, bam_hdr_t* header, bool is_vcf_model);
pVariant str2variant(const std::string& s, pMemoryPool pool);

void format_double(double d, const char* format, std::string& s);

/*!
 * @brief qual 输出格式
 */
std::string& format_qual_value(double d, std::string& s);

void format_vcf_double(double d, std::string& s);
std::string& write_int_vector(const Int32Vector& data, std::string& s);
std::string& write_long_vector(const LongVector& data, std::string& s);
std::string& write_double_vector(const DoubleVector& data, std::string& s);

std::pmr::map<pAllele, char> build_allele_map(const AlleleVector& alleles);
std::pmr::map<pAllele, int32_t> bcf_build_allele_map(const AlleleVector& alleles);

const char* chr_id2name(int32_t id, bam_hdr_t* header);
int32_t chr_name2id(const char* name, bam_hdr_t* header);

void print_reads(const ReadVector& rs, bam_hdr_t* hdr, pMemoryPool pool);
void print_reads(const ReadHashSet& rs, bam_hdr_t* hdr, pMemoryPool pool);
void print_reads(const std::vector<bam1_t*>& rs, bam_hdr_t* hdr, pMemoryPool pool);
void print_reads_to_file(const ReadVector& rs, bam_hdr_t* hdr, pMemoryPool pool);
void print_reads_to_file(const ReadHashSet& rs, bam_hdr_t* hdr, pMemoryPool pool);
void print_reads_to_file(const std::vector<bam1_t*>& rs, bam_hdr_t* hdr, pMemoryPool pool);

void print_haplotypes(const HaplotypeVector& hs, bam_hdr_t* hdr);
void print_haplotypes_to_file(const HaplotypeVector& hs, bam_hdr_t* hdr);

}  // namespace DebugUtils

}  // namespace rovaca

#endif  // ROVACA_HC_DEBUG_UTILS_H_
