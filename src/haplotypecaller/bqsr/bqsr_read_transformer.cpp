#include "bqsr_read_transformer.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include "apply_bqsr_utils.h"
#include "rovaca_logger.h"

static int create_mask(int context_size)
{
    int mask = 0;
    // create 2*contextSize worth of bits
    for (int i = 0; i < context_size; i++) {
        mask = (mask << 2) | 3;
    }
    // shift 4 bits to mask out the bits used to encode the length
    return mask << LENGTH_BITS;
}

static std::string k_standard_covariate_list = "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate";
static std::map<std::string, bool> k_binary_table = {{"true", true}, {"false", false}};
static const uint32_t n_bqsr_resource = 5;

#ifdef __COUNT_INDEL__
static std::map<std::string, EventType> k_event_table = {{"M", BQSR_MATCH}, {"I", BQSR_INSERT}, {"D", BQSR_DELETION}};
#else
static std::map<std::string, EventType> k_event_table = {{"M", BQSR_MATCH}};
#endif

void BQSRReadTransformer::initialize_argument_table(const BQSRReport& report)
{
    recal_arguments_ = std::make_unique<RecalArguments>();
    recal_arguments_->BAQGOP;
    for (int i = 0; i < report.nrows_; ++i) {
        std::string argument = report.underlying_data_[i][0];
        std::string value = report.underlying_data_[i][1];
        if (value == "null") value = "\0";
        if (argument == "covariate" && !value.empty() && (value != k_standard_covariate_list)) {
            throw std::invalid_argument("Non-standard covariates are not supported.");
        }
        else if (argument == "no_standard_covs" && (k_binary_table[value])) {
            throw std::invalid_argument("Non-standard covariates are not supported.");
        }
        else if (argument == "solid_recal_mode" && (recal_arguments_->SOLID_RECAL_MODE != value)) {
            throw std::invalid_argument("Solid is not supported.");
        }
        else if (argument == "solid_nocall_strategy" && (recal_arguments_->SOLID_NOCALL_STRATEGY != value)) {
            throw std::invalid_argument("Solid is not supported.");
        }
        else if (argument == "mismatches_context_size") {
            recal_arguments_->MISMATCHES_CONTEXT_SIZE = std::stoi(value);
        }
        else if (argument == "indels_context_size") {
            recal_arguments_->INDELS_CONTEXT_SIZE = std::stoi(value);
        }
        else if (argument == "mismatches_default_quality") {
            recal_arguments_->MISMATCHES_DEFAULT_QUALITY = std::stoi(value);
        }
        else if (argument == "insertions_default_quality") {
            recal_arguments_->INSERTIONS_DEFAULT_QUALITY = std::stoi(value);
        }
        else if (argument == "deletions_default_quality") {
            recal_arguments_->DELETIONS_DEFAULT_QUALITY = std::stoi(value);
        }
        else if (argument == "maximum_cycle_value") {
            recal_arguments_->MAXIMUM_CYCLE_VALUE = std::stoi(value);
        }
        else if (argument == "low_quality_tail") {
            recal_arguments_->LOW_QUAL_TAIL = std::stoi(value);
        }
        else if (argument == "quantizing_levels") {
            recal_arguments_->QUANTIZING_LEVELS = std::stoi(value);
        }
    }
}

void BQSRReadTransformer::initialize_quantization_table(const BQSRReport& report)
{
    quantization_info_ = std::make_unique<QuantizationInfo>();
    quantization_info_->quantizationLevels_ = recal_arguments_->QUANTIZING_LEVELS;
    quantization_info_->quantized_quals_ = new uint8_t[BQSR_MAX_PHRED_SCORE + 1];
    quantization_info_->empirical_qual_counts_ = new int64_t[BQSR_MAX_PHRED_SCORE + 1];
    const int& qual_col_index = report.col_names_to_index_.at(RecalUtils::QUANTIZED_VALUE_COLUMN_NAME);
    const int& count_col_index = report.col_names_to_index_.at(RecalUtils::QUANTIZED_COUNT_COLUMN_NAME);
    for (int i = 0; i < report.nrows_; ++i) {
        quantization_info_->quantized_quals_[i] = (uint8_t)std::stoi(report.underlying_data_[i][qual_col_index]);
        quantization_info_->empirical_qual_counts_[i] = std::stol(report.underlying_data_[i][count_col_index]);
    }
}

void BQSRReadTransformer::initialize_standard_covariate_list()
{
    covariates_ = std::make_unique<StandardCovariates>();
    for (size_t i = 0; i < read_groups_.size(); ++i) {
        covariates_->rg_covariate_.rg_table_.insert({read_groups_[i], i});
        covariates_->rg_covariate_.rg_reverse_table_.insert({i, read_groups_[i]});
    }

    auto& context_covs = covariates_->context_covariates_;
    context_covs.mismatchesContextSize = recal_arguments_->MISMATCHES_CONTEXT_SIZE;
    context_covs.indelsContextSize = recal_arguments_->INDELS_CONTEXT_SIZE;
    if (context_covs.mismatchesContextSize > CONTEXT_COVS_MAX_DNA_CONTEXT) {
        throw std::invalid_argument("context size cannot be bigger than 13.");
    }
    if (context_covs.indelsContextSize > CONTEXT_COVS_LENGTH_MASK) {
        throw std::invalid_argument("context size cannot be bigger than 15.");
    }
    if (context_covs.mismatchesContextSize <= 0 || context_covs.indelsContextSize <= 0) {
        throw std::invalid_argument("Context size must be positive.");
    }
    context_covs.low_qual_tail = recal_arguments_->LOW_QUAL_TAIL;
    context_covs.mismatchesKeyMask = create_mask(context_covs.mismatchesContextSize);
    context_covs.indelsKeyMask = create_mask(context_covs.indelsContextSize);

    covariates_->cycle_covariates_.MAXIMUM_CYCLE_VALUE = recal_arguments_->MAXIMUM_CYCLE_VALUE;
}

void BQSRReadTransformer::initialize_recalibration_tables(const std::map<std::string, BQSRReport>& tables)
{
    recalibration_tables_ = std::make_unique<RecalibrationTables>();
    int rg_size = static_cast<int>(read_groups_.size());
    int qual_dims = BQSR_MAX_PHRED_SCORE + 1;
    int context_max_key = covariates_->context_covariates_.max_key_value() + 1;
    int cycle_max_key = covariates_->cycle_covariates_.max_key_value() + 1;
    recalibration_tables_->read_group_tables_.create({rg_size, EVENT_TYPE_MAX});
    recalibration_tables_->quality_score_tables_.create({rg_size, qual_dims, EVENT_TYPE_MAX});
    recalibration_tables_->context_tables_.create({rg_size, qual_dims, context_max_key, EVENT_TYPE_MAX});
    recalibration_tables_->cycle_tables_.create({rg_size, qual_dims, cycle_max_key, EVENT_TYPE_MAX});
    const auto& read_group_table = tables.at(RecalUtils::READGROUP_REPORT_TABLE_TITLE);
    const auto& qual_score_table = tables.at(RecalUtils::QUALITY_SCORE_REPORT_TABLE_TITLE);
    const auto& external_table = tables.at(RecalUtils::ALL_COVARIATES_REPORT_TABLE_TITLE);
    parse_read_group_table(recalibration_tables_->read_group_tables_, read_group_table);
    parse_quality_score_table(recalibration_tables_->quality_score_tables_, qual_score_table);
    parse_external_table(external_table);
}

BQSRReadTransformer::BQSRReadTransformer(const std::vector<std::string>& read_groups, const std::string& recal_table_file,
                                         RingMemPool<bam1_t>* mem)
    : read_groups_(read_groups)
    , mem_(mem)
{
    preserve_q_less_than = PRESERVE_QSCORES_LESS_THAN;
    global_qual_score_prior = -1.0;
    emit_original_quals = false;
    use_origin_base_quals = false;
    static_quantized_mapping = nullptr;
    load_report(recal_table_file);
    quantization_info_->non_quantization();

    tasks_resource_ = new BlockingQueue<std::shared_ptr<TransformTask>>(n_bqsr_resource);

    for (uint32_t i = 0; i < n_bqsr_resource; ++i) {
        auto batch_reads = std::make_shared<std::list<bam1_t*>>();
        auto task = std::make_shared<TransformTask>();
        task->batch_reads_ = batch_reads;
        tasks_resource_->push(task);
    }
}

BQSRReadTransformer::~BQSRReadTransformer()
{
    delete[] quantization_info_->quantized_quals_;
    delete[] quantization_info_->empirical_qual_counts_;
    delete tasks_resource_;
}

void BQSRReadTransformer::apply(bam1_t* read)
{
    ReadCovariates target_covariate{read};
    compute_covariates(read, &target_covariate);

    const auto& full_read_key_set = target_covariate.keys;
    int rg_key = full_read_key_set[0][0];

    const RecalDatum& empirical_qual_RG = recalibration_tables_->get_read_group_table().get_value(rg_key, BQSR_MATCH);
    uint8_t* quals = target_covariate.read_base_quality;
    int32_t read_len = target_covariate.read_len;
    double epsilon = global_qual_score_prior > 0.0 ? global_qual_score_prior : empirical_qual_RG.estimated_qual;

    const RecalDatumArray& quality_score_table = recalibration_tables_->get_quality_score_table();
    const RecalDatumArray& context_table = recalibration_tables_->get_context_table();
    const RecalDatumArray& cycle_table = recalibration_tables_->get_cycle_table();
    const uint8_t* quantized_quals = quantization_info_->get_quantized_quals();

    for (int offset = 0; offset < read_len; offset++) {
        if (quals[offset] < preserve_q_less_than) continue;
        const auto& key = full_read_key_set[offset];

        const RecalDatum& empirical_qual_QS = quality_score_table.get_value(key[0], key[1], BQSR_MATCH);
        const RecalDatum& empirical_qual_context = context_table.get_value(key[0], key[1], key[2], BQSR_MATCH);
        const RecalDatum& empirical_qual_cycle = cycle_table.get_value(key[0], key[1], key[3], BQSR_MATCH);
        double recaled_qual = hierarchical_bayesian_quality_estimate(epsilon, empirical_qual_RG, empirical_qual_QS, empirical_qual_context,
                                                                     empirical_qual_cycle);
        // recalibrated quality is bound between 1 and MAX_QUAL.
        uint8_t bound_qual = BOUND_QUAL(FAST_ROUND(recaled_qual), RECALDATUM_MAX_RECALIBRATED_Q_SCORE);
        uint8_t recaled_qual_score = quantized_quals[bound_qual];
        quals[offset] = recaled_qual_score;
        // static_quantized_mapping Not supported currently.
        // quals[offset] = static_quantized_mapping[recaled_qual_score];
    }
}

void BQSRReadTransformer::compute_covariates(bam1_t* read, ReadCovariates* target_covariate)
{
    int low_qual_tail = covariates_->context_covariates_.low_qual_tail;
    int mismatches_context_size = covariates_->context_covariates_.mismatchesContextSize;
    int mismatches_key_mask = covariates_->context_covariates_.mismatchesKeyMask;
    target_covariate->read = read;
    target_covariate->negative_strand = bam_is_rev(read);
    target_covariate->read_base_quality = bam_get_qual(read);
    uint8_t* seq = bam_get_seq(read);
    for (int i = 0; i < read->core.l_qseq; ++i) {
        target_covariate->stranded_clipped_bases[i] = seq_nt16_str[bam_seqi(seq, i)];
    }
    auto& read_len_after_clip = target_covariate->read_len_after_clip;
    auto& stranded_clipped_bases = target_covariate->stranded_clipped_bases;

    read_len_after_clip = get_stranded_clipped_bases(read, stranded_clipped_bases, low_qual_tail);
    context_with(target_covariate->mismatch_keys, read->core.l_qseq, stranded_clipped_bases, mismatches_context_size, mismatches_key_mask);
    target_covariate->record_rg_covariate_value();
    target_covariate->record_qual_score_covariate_value();
    target_covariate->record_context_covariate_value();
    target_covariate->record_cycle_covariate_value();
    return;
}

double BQSRReadTransformer::hierarchical_bayesian_quality_estimate(double epsilon, const RecalDatum& empirical_qual_RG,
                                                                   const RecalDatum& empirical_qual_QS,
                                                                   const RecalDatum& empirical_qual_context,
                                                                   const RecalDatum& empirical_qual_cycle)
{
    double delta_Q = __glibc_likely(empirical_qual_RG.valid) ? empirical_qual_RG.get_empirical_qual(epsilon) - epsilon : 0.0;
    epsilon += delta_Q;
    double delta_Q_reported = __glibc_likely(empirical_qual_QS.valid) ? empirical_qual_QS.get_empirical_qual(epsilon) - (epsilon) : 0.0;
    double delta_Q_covs = 0.0;
    double cond_prior = delta_Q_reported + epsilon;
    if (__glibc_likely(empirical_qual_context.valid)) delta_Q_covs += empirical_qual_context.get_empirical_qual(cond_prior) - cond_prior;
    if (__glibc_likely(empirical_qual_cycle.valid)) delta_Q_covs += empirical_qual_cycle.get_empirical_qual(cond_prior) - cond_prior;
    return cond_prior + delta_Q_covs;
}

static void parse_by_cols(const std::string& line, int ncols, std::vector<std::string>& target)
{
    target.resize(ncols);
    std::istringstream iss(line);
    for (int c = 0; c < ncols; ++c) {
        iss >> target[c];
    }
}

void BQSRReadTransformer::load_report(const std::string& recal_table_file)
{
    std::ifstream recal_table_loader(recal_table_file, std::ios::in);
    std::map<std::string, BQSRReport> all_tables;
    std::string line;
    int ntables;
    char table_name_tmp[128], table_des_tmp[1280];
    std::getline(recal_table_loader, line);  // header
    sscanf(line.c_str(), "%*[^:]:%*[^:]:%d", &ntables);
    for (int i = 0; i < ntables; ++i) {
        int ncols, nrows;
        std::getline(recal_table_loader, line);
        sscanf(line.c_str(), "%*[^:]:%*[^:]:%d:%d", &ncols, &nrows);  // table info
        BQSRReport table(ncols, nrows);

        std::getline(recal_table_loader, line);  // table name and describes.
        sscanf(line.c_str(), "%*[^:]:%*[^:]:%[^:]:%[^:]", table_name_tmp, table_des_tmp);
        table.table_name_ = std::string(table_name_tmp);
        table.table_description_ = std::string(table_des_tmp);
        std::getline(recal_table_loader, line);
        parse_by_cols(line, ncols, table.col_names_);
        for (int idx = 0; idx < ncols; ++idx) table.col_names_to_index_.insert({table.col_names_[idx], idx});

        // fill up table.
        for (int k = 0; k < nrows; ++k) {
            std::getline(recal_table_loader, line);
            parse_by_cols(line, ncols, table.underlying_data_[k]);
        }
        std::getline(recal_table_loader, line);  // empty line marking the end.
        all_tables.insert({table.table_name_, table});
    }

    recal_table_loader.close();
    // Initilize must be in order.
    initialize_argument_table(all_tables.at(RecalUtils::ARGUMENT_REPORT_TABLE_TITLE));
    initialize_quantization_table(all_tables.at(RecalUtils::QUANTIZED_REPORT_TABLE_TITLE));
    initialize_standard_covariate_list();  // all_tables.at(RecalUtils::QUALITY_SCORE_REPORT_TABLE_TITLE)
    initialize_recalibration_tables(all_tables);
    return;
}

void BQSRReadTransformer::parse_read_group_table(RecalDatumArray& array, const BQSRReport& table)
{
    int dim1, dim2;
    const int& rg_col = table.col_names_to_index_.at(RecalUtils::READGROUP_COLUMN_NAME);
    const int& event_col = table.col_names_to_index_.at(RecalUtils::EVENT_TYPE_COLUMN_NAME);
    const int& observations_col = table.col_names_to_index_.at(RecalUtils::NUMBER_OBSERVATIONS_COLUMN_NAME);
    const int& error_col = table.col_names_to_index_.at(RecalUtils::NUMBER_ERRORS_COLUMN_NAME);
    const int& estimated_qual_col = table.col_names_to_index_.at(RecalUtils::ESTIMATED_Q_REPORTED_COLUMN_NAME);
    for (int i = 0; i < table.nrows_; ++i) {
        const std::string& rg = table.underlying_data_[i][rg_col];
        const std::string& event = table.underlying_data_[i][event_col];
        dim1 = covariates_->rg_covariate_.rg_table_.at(rg);
        dim2 = k_event_table.at(event);

        long observation = stol(table.underlying_data_[i][observations_col]);
        double error = stod(table.underlying_data_[i][error_col]);
        double estimated_qual = stod(table.underlying_data_[i][estimated_qual_col]);
        RecalDatum datum{observation, error, estimated_qual};
        array.put_value(datum, {dim1, dim2});
    }
}
void BQSRReadTransformer::parse_quality_score_table(RecalDatumArray& array, const BQSRReport& table)
{
    int dim1, dim2, dim3;
    const int& rg_col = table.col_names_to_index_.at(RecalUtils::READGROUP_COLUMN_NAME);
    const int& qual_col = table.col_names_to_index_.at(RecalUtils::QUALITY_SCORE_COLUMN_NAME);
    const int& event_col = table.col_names_to_index_.at(RecalUtils::EVENT_TYPE_COLUMN_NAME);
    const int& observations_col = table.col_names_to_index_.at(RecalUtils::NUMBER_OBSERVATIONS_COLUMN_NAME);
    const int& error_col = table.col_names_to_index_.at(RecalUtils::NUMBER_ERRORS_COLUMN_NAME);

    for (int i = 0; i < table.nrows_; ++i) {
        const std::string& rg = table.underlying_data_[i][rg_col];
        const std::string& qual = table.underlying_data_[i][qual_col];
        const std::string& event = table.underlying_data_[i][event_col];
        dim1 = covariates_->rg_covariate_.rg_table_.at(rg);
        dim2 = covariates_->quality_score_covariate_.key_from_value(qual);  // TODO: May instead by stoi.
        dim3 = k_event_table.at(event);

        long observation = stol(table.underlying_data_[i][observations_col]);
        double error = convert_to_double(table.underlying_data_[i][error_col], 2);
        double estimated_qual = convert_to_double(table.underlying_data_[i][qual_col], 10);
        RecalDatum datum{observation, error, estimated_qual};
        array.put_value(datum, {dim1, dim2, dim3});
    }
}
void BQSRReadTransformer::parse_external_table(const BQSRReport& table)
{
    int dim1 = 0, dim2 = 0, dim3 = 0, dim4 = 0;
    const int& rg_col = table.col_names_to_index_.at(RecalUtils::READGROUP_COLUMN_NAME);
    const int& qual_col = table.col_names_to_index_.at(RecalUtils::QUALITY_SCORE_COLUMN_NAME);
    const int& cov_name_col = table.col_names_to_index_.at(RecalUtils::COVARIATE_NAME_COLUMN_NAME);
    const int& cov_value_col = table.col_names_to_index_.at(RecalUtils::COVARIATE_VALUE_COLUMN_NAME);
    const int& event_col = table.col_names_to_index_.at(RecalUtils::EVENT_TYPE_COLUMN_NAME);
    const int& observations_col = table.col_names_to_index_.at(RecalUtils::NUMBER_OBSERVATIONS_COLUMN_NAME);
    const int& error_col = table.col_names_to_index_.at(RecalUtils::NUMBER_ERRORS_COLUMN_NAME);

    for (int i = 0; i < table.nrows_; ++i) {
        const std::string& rg = table.underlying_data_[i][rg_col];
        const std::string& qual = table.underlying_data_[i][qual_col];
        const std::string& cov_name = table.underlying_data_[i][cov_name_col];
        const std::string& cov_value = table.underlying_data_[i][cov_value_col];
        const std::string& event = table.underlying_data_[i][event_col];

        dim1 = covariates_->rg_covariate_.rg_table_.at(rg);
        dim2 = covariates_->quality_score_covariate_.key_from_value(qual);  // TODO: May instead by stoi.
        if (cov_name == "Cycle") dim3 = covariates_->cycle_covariates_.key_from_value(cov_value);
        if (cov_name == "Context") dim3 = covariates_->context_covariates_.key_from_value(cov_value);
        dim4 = k_event_table.at(event);

        long observation = stol(table.underlying_data_[i][observations_col]);
        double error = convert_to_double(table.underlying_data_[i][error_col], 2);
        double estimated_qual = convert_to_double(table.underlying_data_[i][qual_col], 10);
        RecalDatum datum{observation, error, estimated_qual};
        recalibration_tables_->get_table_by_name(cov_name).put_value(datum, {dim1, dim2, dim3, dim4});
    }
}

static const uint32_t BQSR_BATCH_READS = (1024 * 1024);
uint32_t reads_in_mem = 0;

bam1_t* bam_dup2(const bam1_t* bsrc, RingMemPool<bam1_t>* mem)
{
    bam1_t* bdst = mem->alloc();

    if (__glibc_unlikely(bdst == nullptr)) {
        bdst = bam_init1();
        bdst = bam_copy1(bdst, bsrc);
        return bdst;
    }
    reads_in_mem++;
    memcpy(bdst->data, bsrc->data, bsrc->l_data);          // copy var-len data
    memcpy(&bdst->core, &bsrc->core, sizeof(bsrc->core));  // copy the rest

    bdst->l_data = bsrc->l_data;
    bdst->m_data = bsrc->m_data;
    bdst->id = bsrc->id;
    return bdst;
}

void BQSRReadTransformer::launch(ReadStream* stream, boost::asio::thread_pool* thread_pool,
                                 BlockingQueue<std::shared_ptr<TransformTask>>* task_queue)
{
    std::atomic_bool reads_finished_{false};
    std::atomic_int tasks_in_progress{0};
    std::mutex mtx;
    std::condition_variable cv;
    uint32_t task_id = 0;
    int32_t last_tid = INT32_MAX;
    long total_count = 0;
    int file_index;

    while (!reads_finished_.load(std::memory_order_acquire)) {
        auto task = tasks_resource_->pop();
        auto& batch_reads = task->batch_reads_;
        task->task_id = task_id;
        uint32_t count = 0;

        do {
            if (__glibc_unlikely(!stream->has_next())) {
                reads_finished_.store(true, std::memory_order_release);
                break;
            }
            bam1_t* item = stream->next(file_index);
            if (item == NULL || item->data == NULL) continue;

            bam1_t* copy = bam_dup2(item, mem_);
            stream->mem_recovery(item);
            batch_reads->emplace_back(copy);
            count++;
            total_count++;

            if (__glibc_unlikely(last_tid != item->core.tid && last_tid != INT32_MAX)) {
                last_tid = item->core.tid;
                break;
            }

        } while (count < BQSR_BATCH_READS);

        tasks_in_progress.fetch_add(1, std::memory_order_relaxed);

        boost::asio::post(*thread_pool, [=, &tasks_in_progress, &mtx, &cv]() mutable {
            pthread_setname_np(pthread_self(), "Covariate");
            std::for_each(batch_reads->begin(), batch_reads->end(), [this](bam1_t* item) { apply(item); });

            task->batch_reads_ = batch_reads;
            task->task_id = task_id;
            task_queue->push(task);

            tasks_in_progress.fetch_sub(1, std::memory_order_relaxed);
            {
                std::unique_lock<std::mutex> lock(mtx);
                cv.notify_all();
            }
        });

        task_id++;
    }
    {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&tasks_in_progress]() { return tasks_in_progress.load(std::memory_order_relaxed) == 0; });
    }
    reads_done_.store(true, std::memory_order_release);
    RovacaLogger::info("ApplyBQSR task count {}", task_id);
    RovacaLogger::info("ApplyBQSR reads count {}", total_count);
    RovacaLogger::info("reads in mem count {}", reads_in_mem);
}

void BQSRReadTransformer::reduce(BlockingQueue<std::shared_ptr<TransformTask>>* task_queue, ReducedQueue* reduced_queue)
{
    long read_count = 0;
    uint32_t current_id{0};
    std::shared_ptr<TransformTask> task_item;
    std::unordered_map<uint32_t, std::shared_ptr<TransformTask>> cache;
    std::unordered_map<uint32_t, std::shared_ptr<TransformTask>>::iterator it;

    while (true) {
        task_item = task_queue->pop(1);

        if (task_item) {
            cache.insert({task_item->task_id, task_item});
            it = cache.find(current_id);

            while (it != cache.end()) {
                auto task = it->second;
                auto& ret = task->batch_reads_;

                cache.erase(it);
                read_count += ret->size();
                reduced_queue->splice(*ret.get());

                tasks_resource_->push(task);

                current_id++;
                it = cache.find(current_id);
            }
        }
        else if (reads_done_.load(std::memory_order_acquire)) {
            reduced_queue->emit_done_signal();
            break;
        }
    }
    RovacaLogger::info("ApplyBQSR done, all reduced reads count {}", read_count);
}
