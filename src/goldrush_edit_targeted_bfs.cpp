#include "utils.hpp"
#include "seqindex.hpp"
#include "mappings.hpp"

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
static const std::string BATCH_MAPPED_IDS_INPUT_READY_PIPE = "batch_mapped_ids_input_ready";
static const std::string MAPPED_IDS_INPUT_PIPE = "mapped_ids_input";
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
serve_batch(const std::string& batch_name,
            const std::string& batch_id_pipe_path,
            const std::string& batch_confirm_pipe_path,
            const std::vector<std::string>& bf_paths_base,
            const std::vector<unsigned>& k_values,
            const size_t cbf_bytes,
            const size_t bf_bytes,
            const unsigned kmer_threshold,
            const unsigned hash_num,
            const std::string& target_seqs_filepath,
            const Index& target_seqs_index,
            const AllMappings& all_mappings,
            const std::string& mapped_seqs_filepath,
            const Index& mapped_seqs_index)
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
  std::vector<std::string> bf_paths;
  for (const auto& bf_path_suffix : bf_paths_base) {
    bf_paths.push_back(batch_prefix + SEPARATOR + bf_path_suffix);
  }

  std::string target_seq_id;
  std::ifstream inputstream(batch_inpipe_path);
  while (inputstream >> contig_id && contig_id != END_SYMBOL) {
    const auto [contig_seq, contig_len] =
      get_seq_with_index<1>(contig_id, contigs_index, contigs_filepath);
    (void)contig_seq;

    const auto contigs_reads_it = contigs_reads.find(contig_id);
    if ((contigs_reads_it == contigs_reads.end()) ||
        (contigs_reads_it->second.empty())) {
      continue;
    }

    const auto& contig_reads_vector = contigs_reads_it->second;
    const auto contig_reads_num = contig_reads_vector.size();
    const auto contig_reads_num_adjusted =
      std::min(contig_reads_num,
               decltype(contig_reads_num)(
                 (double(contig_len) * SUBSAMPLE_MAX_READS_PER_CONTIG_10KBP) /
                 10'000.0));

    const auto random_indices =
      get_random_indices(contig_reads_num, contig_reads_num_adjusted);

    for (const auto read_id_idx : random_indices)
#pragma omp task firstprivate(read_id_idx)                                     \
  shared(bfs, cbfs, contig_reads_vector, reads_index, reads_filepath, ks)
    {
      const auto read_id = contig_reads_vector[read_id_idx].seq_id;
      const auto [seq, seq_len] =
        get_seq_with_index<2>(read_id, reads_index, reads_filepath);
      fill_bfs(seq, seq_len, hash_num, ks, kmer_threshold, cbfs, bfs);
    }
#pragma omp taskwait
  }
  inputstream.close();

  for (size_t i = 0; i < bfs.size(); i++) {
    bfs[i]->save(bf_paths[i]);
  }

  confirm_pipe(batch_confirm_pipepath);

  std::remove(batch_input_pipepath.c_str());
  std::remove(batch_confirm_pipepath.c_str());
}

bool
process_batch_name(const SeqIndex& target_seqs_index,
      const SeqIndex& mapped_seqs_index,
      const AllMappings& all_mappings,
      const size_t cbf_bytes,
      const size_t bf_bytes,
      const unsigned kmer_threshold,
      const std::string& batch_name_input_pipe,
      const std::string& batch_mapped_ids_input_ready_pipe,
      const std::string& mapped_ids_input_pipe,
      const std::string& bfs_ready_pipe,
      const unsigned hash_num,
      const std::vector<unsigned>& k_values) {
  std::string batch_name;
  std::ifstream batch_name_input(batch_name_input_pipe);

  if (!(batch_name_input >> batch_name)) {
    return false;
  }
  batch_name_input.close();
  if (batch_name == END_SYMBOL) {
    return false;
  }

    const auto batch_mapped_ids_input_pipe = batch_name + SEPARATOR +
                                      mapped_ids_input_pipe;
    const auto batch_bfs_ready_pipe = batch_name + SEPARATOR + bfs_ready_pipe;

    make_pipe(batch_mapped_ids_input_pipe);
    make_pipe(batch_bfs_ready_pipe);

    confirm_pipe(batch_mapped_ids_input_ready_pipe);

#pragma omp task firstprivate(batch_name,                                    \
                              batch_mapped_ids_input_pipe,                            \
                              batch_mapped_ids_input_ready_pipe) shared(bf_names,    \
                                                             k_values,               \
                                                             target_seqs_index,    \
                                                             all_mappings,    \
                                                             mapped_seqs_index)
    serve_batch(batch_name,
                batch_mapped_ids_input_pipe,
                batch_mapped_ids_input_ready_pipe,
                bf_names,
                k_values,
                cbf_bytes,
                bf_bytes,
                kmer_threshold,
                hash_num,
                target_seqs_filepath,
                target_seqs_index,
                all_mappings,
                mapped_seqs_filepath,
                mapped_seqs_index);

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
      const std::string& batch_mapped_ids_input_ready_pipe,
      const std::string& mapped_ids_input_pipe,
      const std::string& bfs_ready_pipe,
      const unsigned hash_num,
      const std::vector<unsigned>& k_values)
{
  btllib::check_error(kmer_threshold == 0, FN_NAME + ": k-mer threshold must be >0.");

  make_pipe(batch_name_input_pipe);
  make_pipe(batch_mapped_ids_input_ready_pipe);

  std::vector<std::string> bf_names;
  for (const auto k : ks) {
    bf_names.push_back("k" + std::to_string(k) + BF_EXTENSION);
  }

  btllib::log_info(FN_NAME + ": Accepting batch names at " +
                   batch_name_input_pipepath);
#pragma omp parallel
#pragma omp single
  while (process_batch_name(target_seqs_index,
      mapped_seqs_index,
      all_mappings,
      cbf_bytes,
      bf_bytes,
      kmer_threshold,
      batch_name_input_pipe,
      batch_mapped_ids_input_ready_pipe,
      mapped_ids_input_pipe,
      bfs_ready_pipe,
      hash_num,
      k_values)) {}

  std::remove(batch_name_input_pipe.c_str());
  std::remove(batch_mapped_ids_input_ready_pipe.c_str());

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
        BATCH_MAPPED_IDS_INPUT_READY_PIPE,
        MAPPED_IDS_INPUT_PIPE,
        BFS_READY_PIPE,
        opt::hash_num,
        opt::k_values);

  return 0;
}