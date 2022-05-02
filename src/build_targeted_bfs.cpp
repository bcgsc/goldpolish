#include "utils.hpp"

#include "btllib/status.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/util.hpp"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <memory>
#include <thread>
#include <chrono>
#include <utility>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

static const double SUBSAMPLE_MAX_READS_PER_CONTIG_10KBP = 120.0;
static const double MX_MAX_READS_PER_CONTIG_10KBP = 150.0;
static const unsigned MX_THRESHOLD_MIN = 1;
static const unsigned MX_THRESHOLD_MAX = 30;
static const std::string INPIPE = "input";
static const std::string CONFIRMPIPE = "confirm";
static const std::string SEPARATOR = "_";
static const std::string BF_EXTENSION = ".bf";
static const std::string END_SYMBOL = "x";

namespace opt {
  std::string contigs_filepath;
  std::string contigs_index_filepath;
  std::string contigs_reads_filepath;
  std::string reads_filepath;
  std::string reads_index_filepath;
  size_t cbf_bytes = 10ULL * 1024ULL * 1024ULL;
  size_t bf_bytes = 512ULL * 1024ULL;
  unsigned kmer_threshold = 6;
  std::string prefix = "targeted";
  std::vector<unsigned> ks = { 32, 28, 24, 20 };
  bool ks_set = false;
  unsigned hash_num = 4;
  unsigned threads = 7;
}

void
serve_set(const std::string& set_prefix,
          const std::string& set_input_pipepath,
          const std::string& set_confirm_pipepath,
          const std::vector<std::string>& bf_paths_base,
          const std::vector<unsigned>& ks,
          const size_t cbf_bytes,
          const size_t bf_bytes,
          const unsigned kmer_threshold,
          const unsigned hash_num,
          const std::string& contigs_filepath,
          const Index& contigs_index,
          const ReadsMapping& contigs_reads,
          const std::string& reads_filepath,
          const Index& reads_index)
{
  std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>> cbfs;
  std::vector<std::unique_ptr<btllib::KmerBloomFilter>> bfs;

  for (const auto k : ks) {
    cbfs.push_back(std::make_unique<btllib::KmerCountingBloomFilter8>(cbf_bytes, hash_num, k));
    bfs.push_back(std::make_unique<btllib::KmerBloomFilter>(bf_bytes, hash_num, k));
  }

  std::vector<std::string> bf_paths;
  for (const auto& bf_path_suffix : bf_paths_base) {
    bf_paths.push_back(set_prefix + SEPARATOR + bf_path_suffix);
  }

  std::string contig_id;
  std::ifstream inputstream(set_input_pipepath);
  while (inputstream >> contig_id && contig_id != END_SYMBOL) {
    const auto [contig_seq, contig_len] = get_seq_with_index<1>(contig_id, contigs_index, contigs_filepath);

    const auto contigs_reads_it = contigs_reads.find(contig_id);
    if (
      (contigs_reads_it == contigs_reads.end()) ||
      (contigs_reads_it->second.empty())
        )
    {
      continue;
    }

    const auto& contig_reads_vector = contigs_reads_it->second;
    const auto contig_reads_num = contig_reads_vector.size();
    const auto contig_reads_num_adjusted = std::min(contig_reads_num, decltype(contig_reads_num)((double(contig_len) * SUBSAMPLE_MAX_READS_PER_CONTIG_10KBP) / 10'000.0));

    const auto random_indices = get_random_indices(contig_reads_num, contig_reads_num_adjusted);

    for (const auto read_id_idx : random_indices)
#pragma omp task firstprivate(read_id_idx) shared(bfs, cbfs, contig_reads_vector, reads_index, reads_filepath, ks)
    {
      const auto read_id = contig_reads_vector[read_id_idx].seq_id;
      const auto [seq, seq_len] = get_seq_with_index<2>(read_id, reads_index, reads_filepath);
      fill_bfs(seq, seq_len, hash_num, ks, kmer_threshold, cbfs, bfs);
    }
#pragma omp taskwait
  }
  inputstream.close();

  for (size_t i = 0; i < bfs.size(); i++) {
    bfs[i]->save(bf_paths[i]);
  }

  std::ofstream outputstream(set_confirm_pipepath);
  outputstream << "1" << std::endl;
  outputstream.close();

  std::remove(set_input_pipepath.c_str());
  std::remove(set_confirm_pipepath.c_str());
}

void serve(const std::string& contigs_filepath,
           const Index& contigs_index,
           const ReadsMapping& contigs_reads,
           const std::string& reads_filepath,
           const Index& reads_index,
           const std::string& bf_prefix,
           const size_t cbf_bytes,
           const size_t bf_bytes,
           const unsigned kmer_threshold,
           const std::string& input_pipepath,
           const std::string& confirm_pipepath,
           const unsigned hash_num,
           const std::vector<unsigned>& ks) {
  btllib::check_error(kmer_threshold == 0, "k-mer threshold must be >0.");
  btllib::check_error(mkfifo(input_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");
  btllib::check_error(mkfifo(confirm_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");

  std::vector<std::string> bf_paths_base;
  for (const auto k : ks) {
    bf_paths_base.push_back(std::string(bf_prefix) + "k" + std::to_string(k) + BF_EXTENSION);
  }

  const auto input_pipepath_dirname = btllib::get_dirname(input_pipepath);
  const auto input_pipepath_basename = btllib::get_basename(input_pipepath);
  const auto confirm_pipepath_dirname = btllib::get_dirname(confirm_pipepath);
  const auto confirm_pipepath_basename = btllib::get_basename(confirm_pipepath);

  std::string set_prefix;

  btllib::log_info(std::string("Accepting contig requests at ") + input_pipepath);
#pragma omp parallel
#pragma omp single
  while (true) {
    std::ifstream inputstream(input_pipepath);

    if (!(inputstream >> set_prefix)) { break; }
    inputstream.close();
    if (set_prefix == END_SYMBOL) { break; }

    const auto set_input_pipepath = input_pipepath_dirname + '/' + set_prefix + SEPARATOR + input_pipepath_basename;
    const auto set_confirm_pipepath = confirm_pipepath_dirname + '/' + set_prefix + SEPARATOR + confirm_pipepath_basename;

    btllib::check_error(mkfifo(set_input_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");
    btllib::check_error(mkfifo(set_confirm_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");

    std::ofstream outputstream(confirm_pipepath);
    outputstream << "1" << std::endl;
    outputstream.close();

#pragma omp task firstprivate(set_prefix, set_input_pipepath, set_confirm_pipepath) shared(bf_paths_base, ks, contigs_filepath, contigs_index, contigs_reads, reads_filepath, reads_index)
    serve_set(set_prefix, set_input_pipepath, set_confirm_pipepath, bf_paths_base, ks, cbf_bytes, bf_bytes, kmer_threshold, hash_num, contigs_filepath, contigs_index, contigs_reads, reads_filepath, reads_index);
  }
  btllib::log_info("Targeted BF builder done!");

  std::remove(input_pipepath.c_str());
  std::remove(confirm_pipepath.c_str());
}

int main(int argc, char** argv) {
  btllib::check_error(argc != 7, "Wrong args.");

  bind_to_parent();

  unsigned arg = 1;
  opt::contigs_filepath = argv[arg++];
  opt::contigs_index_filepath = argv[arg++];
  opt::contigs_reads_filepath = argv[arg++];
  opt::reads_filepath = argv[arg++];
  opt::reads_index_filepath = argv[arg++];
  opt::kmer_threshold = std::stoi(argv[arg++]);

  omp_set_nested(1);
  omp_set_num_threads(opt::threads);

  Index reads_index, contigs_index;
  load_index(contigs_index, opt::contigs_index_filepath);
  load_index(reads_index, opt::reads_index_filepath);

  ReadsMapping contigs_reads;
  if (btllib::endswith(opt::contigs_reads_filepath, ".sam") || btllib::endswith(opt::contigs_reads_filepath, ".bam")) {
    load_reads_mapping_sam(contigs_reads, opt::contigs_reads_filepath, contigs_index);
  } else {
    load_reads_mapping(contigs_reads, opt::contigs_reads_filepath, contigs_index, MX_THRESHOLD_MIN);
    filter_read_mappings(contigs_reads, MX_MAX_READS_PER_CONTIG_10KBP, MX_THRESHOLD_MIN, MX_THRESHOLD_MAX, contigs_index);
  }

  serve(opt::contigs_filepath,
        contigs_index,
        contigs_reads,
        opt::reads_filepath,
        reads_index,
        opt::prefix + SEPARATOR,
        opt::cbf_bytes,
        opt::bf_bytes,
        opt::kmer_threshold,
        opt::prefix + SEPARATOR + INPIPE,
        opt::prefix + SEPARATOR + CONFIRMPIPE,
        opt::hash_num,
        opt::ks);

  return 0;
}