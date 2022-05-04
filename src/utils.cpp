#include "utils.hpp"

#include "btllib/status.hpp"
#include "btllib/nthash.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/data_stream.hpp"

#include <algorithm>
#include <cstdlib>
#include <random>
#include <fstream>
#include <thread>
#include <chrono>
#include <memory>
#include <vector>
#include <utility>
#include <string>
#include <limits>
#include <sstream>
#include <cstdio>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static const pid_t INIT_PID = 1;
static const unsigned PARENT_QUERY_PERIOD = 1; // seconds
static const auto FIFO_FLAGS = S_IRUSR | S_IWUSR;

std::vector<size_t> get_random_indices(const size_t total_size, const size_t count) {
  btllib::check_error(count > total_size, FN_NAME + ": count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);
  decltype(all_indices) indices(all_indices.begin(), all_indices.begin() + count);

  return indices;
}

void wait_till_parent_ends() {
  while (getppid() != INIT_PID) {
    std::this_thread::sleep_for(std::chrono::seconds(PARENT_QUERY_PERIOD));
  }
}

void bind_to_parent() {
  (new std::thread([] () {
    wait_till_parent_ends();
    std::exit(EXIT_FAILURE);
  }))->detach();
}

void make_pipe(const std::string& pipepath) {
  btllib::check_error(
      mkfifo(pipepath.c_str(), FIFO_FLAGS) != 0,
      FN_NAME + ": mkfifo failed.");
}

void confirm_pipe(const std::string& pipepath) {
  std::ofstream confirm(pipepath);
  confirm << "1" << std::endl;
}

void fill_bfs(const char* seq,
              const size_t seq_len,
              const unsigned hash_num,
              const std::vector<unsigned>& k_values,
              const unsigned kmer_threshold,
              std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
              std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  for (size_t i = 0; i < k_values.size(); i++)
  {
    const auto k = k_values[i];
    const auto& cbf = cbfs[i];
    const auto& bf = bfs[i];
    btllib::NtHash nthash(seq, seq_len, hash_num, k);
    while (nthash.roll()) {
      if (cbf->insert_thresh_contains(nthash.hashes(), kmer_threshold) >= kmer_threshold) {
        bf->insert(nthash.hashes());
      }
    }
  }
}