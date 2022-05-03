#include "common.hpp"

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
  btllib::check_error(count > total_size, FN_NAME + ": count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);
  decltype(all_indices) indices(all_indices.begin(), all_indices.begin() + count);

  return indices;
}

void load_index(Index& index, const std::string& filepath) {
  btllib::log_info(FN_NAME + ": Loading index from " + filepath + "... ");
  
  // Parameter sanity check
  btllib::check_error(!index.empty(), FN_NAME + ": index is not empty.");

  std::ifstream ifs(filepath);
  std::string token, id;
  unsigned long /*id_start = 0,*/ id_end = 0, seq_end = 0;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 4) {
      case 0: id = std::move(token); break;
      case 1: /*id_start = std::stoul(token);*/ break;
      case 2: id_end = std::stoul(token); break;
      case 3: {
        seq_end = std::stoul(token);
        index.emplace(std::piecewise_construct, std::make_tuple(id), std::make_tuple(id_end + 1, seq_end - id_end - 1));
        break;
      }
      default: {
        btllib::log_error(FN_NAME + ": Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info(FN_NAME + ": Done!");
}

void
load_mappings(AllMappings& all_mappings, const std::string& filepath, const Index& target_seqs_index, const unsigned mx_threshold_min, const unsigned mx_threshold_max, const double mx_max_mapped_seqs_per_target_10kbp) {
  if (btllib::endswith(opt::filepath, ".sam") || btllib::endswith(opt::filepath, ".bam")) {
    load_mappings_sam(mappings, opt::filepath, targed_seqs_index);
  } else {
    load_mappings_ntlink(mappings, opt::filepath, target_seqs_index, mx_threshold_min);
    filter_mappings(mappings, mx_max_mapped_seqs_per_target_10kbp, mx_threshold_min, mx_threshold_max, target_seqs_index);
  }
}

void load_mappings_ntlink(AllMappings& all_mappings, const std::string& filepath, const Index& target_seqs_index, const unsigned mx_threshold_min) {
  btllib::log_info(FN_NAME + ": Loading mappings from " + filepath + "... ");

  // Parameter sanity check
  btllib::check_error(!all_mappings.empty(), FN_NAME + ": all_mappings is not empty.");

  std::ifstream ifs(filepath);
  std::string token, mapped_seq_id, target_seq_id;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 3) {
      case 0: mapped_seq_id = std::move(token); break;
      case 1: target_seq_id = std::move(token); break;
      case 2: {
        if (target_seqs_index.find(target_seq_id) == target_seqs_index.end()) {
          break;
        }
        auto it = all_mappings.find(target_seq_id);
        if (it == all_mappings.end()) {
          const auto emplacement = all_mappings.emplace(target_seq_id, Mappings());
          it = emplacement.first;
        }
        const auto minimizers = std::stoul(token);
        if (minimizers >= mx_threshold_min) {
          it->second.push_back(Mapping(mapped_seq_id, minimizers));
        }
        break;
      }
      default: {
        btllib::log_error(FN_NAME + ": Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info(FN_NAME + ": Done!");
}

void
load_mappings_sam(AllMappings& all_mappings, const std::string& filepath, const Index& target_seqs_index) {
  btllib::log_info(FN_NAME + ": Loading mappings from " + filepath + "... ");

  // Parameter sanity check
  btllib::check_error(!mappings.empty(), FN_NAME + ": all_mappings is not empty.");

  size_t n = 1024;
  char* line = new char[n];
  btllib::DataSource data_source(filepath);

  std::string token, mapped_seq_id, target_seq_id;
  unsigned long i = 0;
  while (getline(&line, &n, data_source)) {
    if (n > 0 && line[0] == '@') {
      continue;
    }
    std::stringstream ss(line);
    while (ss >> token) {
      switch (i % 3) {
        case 0: mapped_seq_id = std::move(token); break;
        case 2: target_seq_id = std::move(token); break;
        case 4: {
          if (target_seqs_index.find(target_seq_id) == target_seqs_index.end()) {
            break;
          }
          auto it = mappings.find(contig_id);
          if (it == mappings.end()) {
            const auto emplacement = mappings.emplace(contig_id, Mappings());
            it = emplacement.first;
          }
          it->second.push_back(Mapping(read_id, 0));
          break;
        }
        default: {
          btllib::log_error(FN_NAME + ": Invalid switch branch.");
          std::exit(EXIT_FAILURE);
        }
      }
      i++;
    }
  }
  btllib::log_info(FN_NAME + ": Done!");
}

unsigned mapped_seqs_for_threshold(const mappings& mappings, const unsigned mx_threshold) {
  unsigned mapped_seq_num = 0;
  std::for_each(mappings.begin(), mappings.end(), [&](const Mapping& mapping) {
    if (mapping.mx_in_common >= mx_threshold) {
      mapped_seq_num++;
    }
  });
  return mapped_seq_num;
}

void
filter_mappings(Mappings& mappings, const double max_mapped_seqs_per_target_10kbp, const unsigned mx_threshold_min, const unsigned mx_threshold_max, Index& target_seqs_index) {
  btllib::log_info(FN_NAME + ": Filtering contig mapped_seqs... "));

  // Parameter sanity check
  btllib::check_error(mappings.empty(), FN_NAME + ": mappings is empty.");
  btllib::check_error(target_seqs_index.empty(), FN_NAME + ": target_seqs_index is empty.");
  btllib::check_error(max_mapped_seqs_per_target_10kbp <= 0, FN_NAME + ": max_mapped_seqs_per_target_10kbp is not positive.");
  btllib::check_error(mx_threshold_min >= mx_threshold_max, FN_NAME + ": mx_threshold_min is not smaller than mx_threshold_max.");

  for (auto& target_mappings : mappings) {
    const auto& target_seq_id = target_mappings.first;
    const auto& m = target_mappings.second;
    if (m.empty()) { continue; }

    if (taget_seqs_index.find(target_seq_id) == target_seqs_index.end()) {
      continue;
    }
    const auto target_seq_len = target_seqs_index.at(target_seq_id).seq_len;
    const int max_mapped_seqs = std::ceil(double(target_seq_len) * max_mapped_seqs_per_target_10kbp / 10'000.0);
    btllib::check_error(max_mapped_seqs <= 0, FN_NAME + ": max_mapped_seqs <= 0.");

    int updated_mx_threshold_min = mx_threshold_min;
    int updated_mx_threshold_min_mapped_seqs = m.size();

    int updated_mx_threshold_max = mx_threshold_max;
    int updated_mx_threshold_max_mapped_seqs = mapped_seqs_for_threshold(m, updated_mx_threshold_max);

    int mx_threshold = -1;
    if (updated_mx_threshold_min_mapped_seqs <= max_mapped_seqs) {
      mx_threshold = updated_mx_threshold_min;
    } else if (updated_mx_threshold_max_mapped_seqs > max_mapped_seqs) {
      mx_threshold = updated_mx_threshold_max;
    } else {
      while (updated_mx_threshold_max - updated_mx_threshold_min > 1) {
        const int mx_threshold_mid = (updated_mx_threshold_max + updated_mx_threshold_min) / 2;
        const int mx_threshold_mid_mapped_seqs = mapped_seqs_for_threshold(m, mx_threshold_mid);
        if (mx_threshold_mid_mapped_seqs > max_mapped_seqs) {
          updated_mx_threshold_min = mx_threshold_mid;
          updated_mx_threshold_min_mapped_seqs = mx_threshold_mid_mapped_seqs;
        } else {
          updated_mx_threshold_max = mx_threshold_mid;
          updated_mx_threshold_max_mapped_seqs = mx_threshold_mid_mapped_seqs;
        }
      }
      mx_threshold = updated_mx_threshold_max;

      btllib::check_error(std::abs(updated_mx_threshold_min - updated_mx_threshold_max) > 1, FN_NAME + ": abs(updated_mx_threshold_min - updated_mx_threshold_max) > 1.");
      btllib::check_error(updated_mx_threshold_min_mapped_seqs <= max_mapped_seqs, FN_NAME + ": Min threshold wrongly calculated.");
      btllib::check_error(updated_mx_threshold_max_mapped_seqs > max_mapped_seqs, FN_NAME + ": Max threshold wrongly calculated.");
      btllib::check_error(updated_mx_threshold_min >= updated_mx_threshold_max, FN_NAME + ": updated_mx_threshold_min >= updated_mx_threshold_max.");
    }

    btllib::check_error(mx_threshold < 0, FN_NAME + ": mx_threshold < 0.");
    btllib::check_error(mx_threshold < int(mx_threshold_min), FN_NAME + ": mx_threshold < mx_threshold_min.");
    btllib::check_error(mx_threshold > int(mx_threshold_max), FN_NAME + ": mx_threshold > mx_threshold_max.");

    decltype(target_mappings.second) new_m;
    std::for_each(m.begin(), m.end(), [&](const Mapping& mapping) {
      if (int(mapping.mx_in_common) >= mx_threshold) {
        new_m.push_back(mapping);
      }
    });
    target_mappings.second = new_m;
  }
  btllib::log_info(FN_NAME + ": Done!");
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