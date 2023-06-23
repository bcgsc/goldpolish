#include "utils.hpp"

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static const pid_t INIT_PID = 1;
static const unsigned PARENT_QUERY_PERIOD = 1; // seconds
static const auto FIFO_FLAGS = S_IRUSR | S_IWUSR;

std::vector<size_t>
get_random_indices(const size_t total_size, const size_t count)
{
  btllib::check_error(count > total_size,
                      FN_NAME + ": count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);

  using advance_type = typename std::iterator_traits<
    decltype(all_indices.begin())>::difference_type;

  decltype(all_indices) indices(
    all_indices.begin(), std::next(all_indices.begin(), advance_type(count)));

  return indices;
}

void
wait_till_parent_ends()
{
  while (getppid() != INIT_PID) {
    std::this_thread::sleep_for(std::chrono::seconds(PARENT_QUERY_PERIOD));
  }
}

void
bind_to_parent()
{
  (new std::thread([]() {
    wait_till_parent_ends();
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }))
    ->detach();
}

void
make_pipe(const std::string& pipepath)
{
  const auto ret = mkfifo(pipepath.c_str(), FIFO_FLAGS);
  btllib::check_error(ret != 0,
                      FN_NAME + ": mkfifo failed: " + btllib::get_strerror());
}

std::string
read_pipe(const std::string& pipepath)
{
  std::string str;
  std::ifstream input(pipepath);
  input >> str;
  return str;
}

void
confirm_pipe(const std::string& pipepath)
{
  std::ofstream confirm(pipepath);
  confirm << "1" << std::endl;
}

void
fill_bfs(const char* seq,
         const size_t seq_len,
         const unsigned hash_num,
         const std::vector<unsigned>& k_values,
         const unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  btllib::check_error(kmer_threshold < 4,
                      FN_NAME + ": kmer_threshold must be "
                                "greater than or equal to 4.");
  unsigned adjusted_kmer_threshold = kmer_threshold - 2;
  for (size_t i = 0; i < k_values.size(); i++) {
    const auto k = k_values[i];
    const auto& cbf = cbfs[i];
    const auto& bf = bfs[i];
    btllib::NtHash nthash(seq, seq_len, hash_num, k);
    while (nthash.roll()) {
      if (cbf->insert_thresh_contains(nthash.hashes(),
                                      adjusted_kmer_threshold) >=
          adjusted_kmer_threshold) {
        bf->insert(nthash.hashes());
      }
    }
    adjusted_kmer_threshold++;
  }
}
