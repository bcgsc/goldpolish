#include "mappings.hpp"
#include "seqindex.hpp"
#include "utils.hpp"

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static const double MX_MAX_MAPPED_SEQS_PER_TARGET_10KBP = 150.0;
static const double SUBSAMPLE_MAX_MAPPED_SEQS_PER_TARGET_10KBP = 120.0;
static const unsigned MX_THRESHOLD_MIN = 1;
static const unsigned MX_THRESHOLD_MAX = 30;
static const std::string BATCH_NAME_INPUT_PIPE = "batch_name_input";
static const std::string BATCH_TARGET_IDS_INPUT_READY_PIPE =
  "batch_target_ids_input_ready";
static const std::string TARGET_IDS_INPUT_PIPE = "target_ids_input";
static const std::string BFS_READY_PIPE = "bfs_ready";
static const std::string SEPARATOR = "-";
static const std::string BF_EXTENSION = ".bf";
static const std::string END_SYMBOL = "x";

namespace opt {
std::string target_seqs_filepath;
std::string target_seqs_index_filepath;
std::string mappings_filepath;
std::string mapped_seqs_filepath;
std::string mapped_seqs_index_filepath;
size_t cbf_bytes = 10ULL * 1024ULL * 1024ULL;
size_t bf_bytes = 512ULL * 1024ULL;
unsigned kmer_threshold = 5;
std::vector<unsigned> k_values = { 32, 28, 24, 20 };
unsigned hash_num = 4;
unsigned threads = 7;
}

void
serve_batch(const SeqIndex& target_seqs_index,
            const SeqIndex& mapped_seqs_index,
            const AllMappings& all_mappings,
            const size_t cbf_bytes,
            const size_t bf_bytes,
            const unsigned kmer_threshold,
            const std::string& batch_name,
            const std::string& target_ids_input_pipe,
            const std::string& bfs_ready_pipe,
            const std::vector<std::string>& bf_names,
            const unsigned hash_num,
            const std::vector<unsigned>& k_values,
            const double subsample_max_mapped_seqs_per_target_10kbp
            )
{
  // Initialize Bloom filters
  std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>> cbfs;
  std::vector<std::unique_ptr<btllib::KmerBloomFilter>> bfs;
  for (const auto k : k_values) {
    cbfs.push_back(std::make_unique<btllib::KmerCountingBloomFilter8>(
      cbf_bytes, hash_num, k));
    bfs.push_back(
      std::make_unique<btllib::KmerBloomFilter>(bf_bytes, hash_num, k));
  }

  // Set Bloom filter paths
  std::vector<std::string> bf_full_names;
  for (const auto& bf_name : bf_names) {
    bf_full_names.push_back(batch_name + SEPARATOR + bf_name);
  }

  std::string target_seq_id;
  std::ifstream inputstream(target_ids_input_pipe);
  while (inputstream >> target_seq_id && target_seq_id != END_SYMBOL) {
    const auto target_seq_len = target_seqs_index.get_seq_len(target_seq_id);

    const auto& mappings = all_mappings.get_mappings(target_seq_id);
    if (mappings.empty()) {
      continue;
    }
    const auto mappings_num = mappings.size();
    const auto mappings_num_adjusted = std::min(mappings_num,
               decltype(mappings_num)(
                 (double(target_seq_len) * subsample_max_mapped_seqs_per_target_10kbp) /
                 10'000.0));

    const auto random_indices =
      get_random_indices(mappings_num, mappings_num_adjusted);

    for (const auto mapped_id_idx : random_indices)
#pragma omp task firstprivate(mapped_id_idx)                                     \
  shared(bfs, cbfs, mappings, mapped_seqs_index, k_values)
    {
      const auto mapped_id = mappings[mapped_id_idx].seq_id;
      const auto [seq, seq_len] = mapped_seqs_index.get_seq<1>(mapped_id);
      fill_bfs(seq, seq_len, hash_num, k_values, kmer_threshold, cbfs, bfs);
    }
#pragma omp taskwait
  }
  inputstream.close();

  for (size_t i = 0; i < bfs.size(); i++) {
    bfs[i]->save(bf_full_names[i]);
  }

  confirm_pipe(bfs_ready_pipe);

  std::remove(target_ids_input_pipe.c_str());
  std::remove(bfs_ready_pipe.c_str());
}

bool
process_batch_name(const SeqIndex& target_seqs_index,
                   const SeqIndex& mapped_seqs_index,
                   const AllMappings& all_mappings,
                   const size_t cbf_bytes,
                   const size_t bf_bytes,
                   const unsigned kmer_threshold,
                   const std::string& batch_name_input_pipe,
                   const std::string& batch_target_ids_input_ready_pipe,
                   const std::string& target_ids_input_pipe,
                   const std::string& bfs_ready_pipe,
                   const unsigned hash_num,
                   const std::vector<unsigned>& k_values,
                   const std::vector<std::string>& bf_names,
                   const double subsample_max_mapped_seqs_per_target_10kbp)
{
  const auto batch_name = read_pipe(batch_name_input_pipe);
  if (batch_name.empty() || batch_name == END_SYMBOL) {
    return false;
  }

  const auto batch_target_ids_input_pipe =
    batch_name + SEPARATOR + target_ids_input_pipe;
  const auto batch_bfs_ready_pipe = batch_name + SEPARATOR + bfs_ready_pipe;

  make_pipe(batch_target_ids_input_pipe);
  make_pipe(batch_bfs_ready_pipe);

  confirm_pipe(batch_target_ids_input_ready_pipe);

#pragma omp task firstprivate(                                                 \
  batch_name, batch_target_ids_input_pipe, batch_target_ids_input_ready_pipe)  \
  shared(                                                                      \
    target_seqs_index, mapped_seqs_index, all_mappings, bf_names, k_values)
  serve_batch(target_seqs_index,
              mapped_seqs_index,
              all_mappings,
              cbf_bytes,
              bf_bytes,
              kmer_threshold,
              batch_name,
              batch_target_ids_input_pipe,
              batch_bfs_ready_pipe,
              bf_names,
              hash_num,
              k_values,
              subsample_max_mapped_seqs_per_target_10kbp);

  return true;
}

void
serve(const SeqIndex& target_seqs_index,
      const SeqIndex& mapped_seqs_index,
      const AllMappings& all_mappings,
      const size_t cbf_bytes,
      const size_t bf_bytes,
      const unsigned kmer_threshold,
      const std::string& batch_name_input_pipe,
      const std::string& batch_target_ids_input_ready_pipe,
      const std::string& target_ids_input_pipe,
      const std::string& bfs_ready_pipe,
      const unsigned hash_num,
      const std::vector<unsigned>& k_values,
      const double subsample_max_mapped_seqs_per_target_10kbp)
{
  btllib::check_error(kmer_threshold == 0,
                      FN_NAME + ": k-mer threshold must be >0.");

  make_pipe(batch_name_input_pipe);
  make_pipe(batch_target_ids_input_ready_pipe);

  std::vector<std::string> bf_names;
  for (const auto k : k_values) {
    bf_names.push_back("k" + std::to_string(k) + BF_EXTENSION);
  }

  btllib::log_info(FN_NAME + ": Accepting batch names at " +
                   batch_name_input_pipe);
#pragma omp parallel
#pragma omp single
  while (process_batch_name(target_seqs_index,
                            mapped_seqs_index,
                            all_mappings,
                            cbf_bytes,
                            bf_bytes,
                            kmer_threshold,
                            batch_name_input_pipe,
                            batch_target_ids_input_ready_pipe,
                            target_ids_input_pipe,
                            bfs_ready_pipe,
                            hash_num,
                            k_values,
                            bf_names,
                            subsample_max_mapped_seqs_per_target_10kbp)) {
  }

  std::remove(batch_name_input_pipe.c_str());
  std::remove(batch_target_ids_input_ready_pipe.c_str());

  btllib::log_info(FN_NAME + ": Targeted BF builder done!");
}

int
main(int argc, char** argv)
{
  btllib::check_error(argc != 7, "Wrong args.");

  bind_to_parent();

  unsigned arg = 1;
  opt::target_seqs_filepath = argv[arg++];
  opt::target_seqs_index_filepath = argv[arg++];
  opt::mappings_filepath = argv[arg++];
  opt::mapped_seqs_filepath = argv[arg++];
  opt::mapped_seqs_index_filepath = argv[arg++];
  opt::kmer_threshold = std::stoi(argv[arg++]);

  omp_set_nested(1);
  omp_set_num_threads(opt::threads);

  SeqIndex target_seqs_index(opt::target_seqs_index_filepath,
                             opt::target_seqs_filepath);
  SeqIndex mapped_seqs_index(opt::mapped_seqs_index_filepath,
                             opt::mapped_seqs_filepath);

  AllMappings all_mappings(opt::mappings_filepath,
                           target_seqs_index,
                           MX_THRESHOLD_MIN,
                           MX_THRESHOLD_MAX,
                           MX_MAX_MAPPED_SEQS_PER_TARGET_10KBP);

  serve(target_seqs_index,
        mapped_seqs_index,
        all_mappings,
        opt::cbf_bytes,
        opt::bf_bytes,
        opt::kmer_threshold,
        BATCH_NAME_INPUT_PIPE,
        BATCH_TARGET_IDS_INPUT_READY_PIPE,
        TARGET_IDS_INPUT_PIPE,
        BFS_READY_PIPE,
        opt::hash_num,
        opt::k_values,
        SUBSAMPLE_MAX_MAPPED_SEQS_PER_TARGET_10KBP);

  return 0;
}