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

const pid_t INIT_PID = 1;
const unsigned PARENT_QUERY_PERIOD = 1; // seconds

std::vector<size_t> get_random_indices(const size_t total_size, const size_t count) {
  btllib::check_error(count > total_size, "get_random_indices: count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);
  decltype(all_indices) indices(all_indices.begin(), all_indices.begin() + count);

  return indices;
}

void load_index(Index& index, const std::string& filepath) {
  btllib::log_info(std::string("load_index: Loading index from ") + filepath + "... ");
  
  // Parameter sanity check
  btllib::check_error(!index.empty(), "load_index: index is not empty.");

  std::ifstream ifs(filepath);
  std::string token, id;
  unsigned long id_start, id_end, seq_end;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 4) {
      case 0: id = std::move(token); break;
      case 1: id_start = std::stoul(token); break;
      case 2: id_end = std::stoul(token); break;
      case 3: {
        seq_end = std::stoul(token);
        index.emplace(std::piecewise_construct, std::make_tuple(id), std::make_tuple(id_end + 1, seq_end - id_end - 1));
        break;
      }
      default: {
        btllib::log_error("load_index: Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("load_index: Done!");
}

void load_reads_mapping(ReadsMapping& reads_mapping, const std::string& filepath, const Index& contigs_index, const unsigned mx_threshold_min) {
  btllib::log_info(std::string("load_reads_mapping: Loading contig reads from ") + filepath + "... ");

  // Parameter sanity check
  btllib::check_error(!reads_mapping.empty(), "load_reads_mapping: reads_mapping is not empty.");

  std::ifstream ifs(filepath);
  std::string token, read_id, contig_id;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 3) {
      case 0: read_id = std::move(token); break;
      case 1: contig_id = std::move(token); break;
      case 2: {
        if (contigs_index.find(contig_id) == contigs_index.end()) {
          break;
        }
        auto it = reads_mapping.find(contig_id);
        if (it == reads_mapping.end()) {
          const auto emplacement = reads_mapping.emplace(contig_id, std::vector<ReadMapping>());
          it = emplacement.first;
        }
        const auto minimizers = std::stoul(token);
        if (minimizers >= mx_threshold_min) {
          it->second.push_back(ReadMapping(read_id, minimizers));
        }
        break;
      }
      default: {
        btllib::log_error("load_reads_mapping: Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("load_reads_mapping: Done!");
}

unsigned reads_for_threshold(const std::vector<ReadMapping>& mappings, const unsigned mx_threshold) {
  unsigned read_num = 0;
  std::for_each(mappings.begin(), mappings.end(), [&](const ReadMapping& mapping) {
    if (mapping.mx_in_common >= mx_threshold) {
      read_num++;
    }
  });
  return read_num;
}

void
filter_read_mappings(ReadsMapping& reads_mapping, const double max_reads_per_contig_10kbp, const unsigned mx_threshold_min, const unsigned mx_threshold_max, Index& contigs_index) {
  btllib::log_info(std::string("filter_read_mappings: Filtering contig reads... "));

  // Parameter sanity check
  btllib::check_error(reads_mapping.empty(), "filter_read_mappings: reads_mapping is empty.");
  btllib::check_error(contigs_index.empty(), "filter_read_mappings: contigs_index is empty.");
  btllib::check_error(max_reads_per_contig_10kbp <= 0, "filter_read_mappings: max_reads_per_contig_10kbp is not positive.");
  btllib::check_error(mx_threshold_min >= mx_threshold_max, "filter_read_mappings: mx_threshold_min is not smaller than mx_threshold_max.");

  for (auto& contig_mappings : reads_mapping) {
    const auto& contig_id = contig_mappings.first;
    const auto& mappings = contig_mappings.second;
    if (mappings.empty()) { continue; }

    if (contigs_index.find(contig_id) == contigs_index.end()) {
      continue;
    }
    const auto contig_len = contigs_index.at(contig_id).seq_len;
    const int max_reads = std::ceil(double(contig_len) * max_reads_per_contig_10kbp / 10'000.0);
    btllib::check_error(max_reads <= 0, "filter_read_mappings: max_reads <= 0.");

    int updated_mx_threshold_min = mx_threshold_min;
    int updated_mx_threshold_min_reads = mappings.size();

    int updated_mx_threshold_max = mx_threshold_max;
    int updated_mx_threshold_max_reads = reads_for_threshold(mappings, updated_mx_threshold_max);

    int mx_threshold = -1;
    if (updated_mx_threshold_min_reads <= max_reads) {
      mx_threshold = updated_mx_threshold_min;
    } else if (updated_mx_threshold_max_reads > max_reads) {
      mx_threshold = updated_mx_threshold_max;
    } else {
      while (updated_mx_threshold_max - updated_mx_threshold_min > 1) {
        const int mx_threshold_mid = (updated_mx_threshold_max + updated_mx_threshold_min) / 2;
        const int mx_threshold_mid_reads = reads_for_threshold(mappings, mx_threshold_mid);
        if (mx_threshold_mid_reads > max_reads) {
          updated_mx_threshold_min = mx_threshold_mid;
          updated_mx_threshold_min_reads = mx_threshold_mid_reads;
        } else {
          updated_mx_threshold_max = mx_threshold_mid;
          updated_mx_threshold_max_reads = mx_threshold_mid_reads;
        }
      }
      mx_threshold = updated_mx_threshold_max;

      btllib::check_error(std::abs(updated_mx_threshold_min - updated_mx_threshold_max) > 1, "filter_read_mappings: abs(updated_mx_threshold_min - updated_mx_threshold_max) > 1.");
      btllib::check_error(updated_mx_threshold_min_reads <= max_reads, "filter_read_mappings: Min threshold wrongly calculated.");
      btllib::check_error(updated_mx_threshold_max_reads > max_reads, "filter_read_mappings: Max threshold wrongly calculated.");
      btllib::check_error(updated_mx_threshold_min >= updated_mx_threshold_max, "filter_read_mappings: updated_mx_threshold_min >= updated_mx_threshold_max.");
    }

    btllib::check_error(mx_threshold < 0, "filter_read_mappings: mx_threshold < 0.");
    btllib::check_error(mx_threshold < mx_threshold_min, "filter_read_mappings: mx_threshold < mx_threshold_min.");
    btllib::check_error(mx_threshold > mx_threshold_max, "filter_read_mappings: mx_threshold > mx_threshold_max.");

    decltype(contig_mappings.second) new_mappings;
    std::for_each(mappings.begin(), mappings.end(), [&](const ReadMapping& mapping) {
      if (mapping.mx_in_common >= mx_threshold) {
        new_mappings.push_back(mapping);
      }
    });
    contig_mappings.second = new_mappings;
  }
  btllib::log_info("filter_read_mappings: Done!");
}

void
load_reads_mapping_sam(ReadsMapping& reads_mapping, const std::string& filepath, const Index& contigs_index) {
  btllib::log_info(std::string("load_reads_mapping_sam: Loading contig reads from ") + filepath + "... ");

  // Parameter sanity check
  btllib::check_error(!reads_mapping.empty(), "load_reads_mapping: reads_mapping is not empty.");

  size_t n = 1024;
  char* line = new char[n];
  btllib::DataSource data_source(filepath);

  std::string token, read_id, contig_id;
  unsigned long i = 0;
  while (getline(&line, &n, data_source)) {
    if (n > 0 && line[0] == '@') {
      continue;
    }
    std::stringstream ss(line);
    while (ss >> token) {
      switch (i % 3) {
        case 0: read_id = std::move(token); break;
        case 2: contig_id = std::move(token); break;
        case 4: {
          if (contigs_index.find(contig_id) == contigs_index.end()) {
            break;
          }
          auto it = reads_mapping.find(contig_id);
          if (it == reads_mapping.end()) {
            const auto emplacement = reads_mapping.emplace(contig_id, std::vector<ReadMapping>());
            it = emplacement.first;
          }
          it->second.push_back(ReadMapping(read_id, 0));
          break;
        }
        default: {
          btllib::log_error("load_reads_mapping_sam: Invalid switch branch.");
          std::exit(EXIT_FAILURE);
        }
      }
      i++;
    }
  }
  btllib::log_info("load_reads_mapping_sam: Done!");
}

void fill_bfs(const char* seq,
              const size_t seq_len,
              const unsigned hash_num,
              const std::vector<unsigned>& ks,
              const unsigned kmer_threshold,
              std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
              std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  for (size_t i = 0; i < ks.size(); i++)
  {
    const auto& cbf = cbfs[i];
    const auto& bf = bfs[i];
    btllib::NtHash nthash(seq, seq_len, hash_num, ks[i]);
    while (nthash.roll()) {
      if (cbf->insert_thresh_contains(nthash.hashes(), kmer_threshold) >= kmer_threshold) {
        bf->insert(nthash.hashes());
      }
    }
  }
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