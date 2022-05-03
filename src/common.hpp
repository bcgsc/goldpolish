#ifndef COMMON_HPP
#define COMMON_HPP

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#define FN_NAME (std::string(__FUNCTION__))

using SeqId = std::string;

struct Mapping
{
  SeqId seq_id;
  unsigned mx_in_common;

  Mapping(const SeqId& seq_id, const unsigned mx_in_common)
    : seq_id(seq_id)
    , mx_in_common(mx_in_common)
  {}
};

struct SeqCoordinates
{
  size_t seq_start, seq_len;

  SeqCoordinates(const size_t seq_start, const size_t seq_len)
    : seq_start(seq_start)
    , seq_len(seq_len)
  {}
};

using Mappings = std::vector<Mapping>;
using AllMappings = std::unordered_map<SeqId, Mappings>;
using Index = std::unordered_map<SeqId, SeqCoordinates>;

std::vector<size_t>
get_random_indices(size_t total_size, size_t count);

void
load_index(Index& index, const std::string& filepath);

void
load_mappings(AllMappings& all_mappings,
              const std::string& filepath,
              const Index& target_seqs_index,
              unsigned mx_threshold_min,
              unsigned mx_threshold_max,
              double mx_max_mapped_seqs_per_target_10kbp);

void
load_mappings_ntlink(AllMappings& all_mappings,
                     const std::string& filepath,
                     const Index& target_seqs_index,
                     unsigned mx_threshold_min);

void
load_mappings_sam(AllMappings& all_mappings,
                  const std::string& filepath,
                  const Index& target_seqs_index);

void
filter_mappings(AllMappings& all_mappings,
                double max_mapped_seqs_per_target_10kbp,
                unsigned mx_threshold_min,
                unsigned mx_threshold_max,
                Index& target_seqs_index);

void
fill_bfs(const char* seq,
         size_t seq_len,
         unsigned hash_num,
         const std::vector<unsigned>& k_values,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs);

inline void
fill_bfs(const std::string& seq,
         unsigned hash_num,
         const std::vector<unsigned>& k_values,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs)
{
  fill_bfs(
    seq.c_str(), seq.size(), hash_num, k_values, kmer_threshold, cbfs, bfs);
}

void
wait_till_parent_ends();

void
bind_to_parent();

template<int i>
std::tuple<const char*, size_t>
get_seq_with_index(const std::string& id,
                   const Index& index,
                   const std::string& seqs_filepath)
{
  thread_local static std::ifstream* seqs_file;
  thread_local static bool seqs_file_initialized = false;

  static const size_t max_seqlen = 1024ULL * 1024ULL;
  thread_local static char* seq;
  thread_local static bool seq_initialized = false;

  if (!seqs_file_initialized) {
    seqs_file = new std::ifstream(seqs_filepath);
    seqs_file_initialized = true;
  }

  if (!seq_initialized) {
    seq = new char[max_seqlen];
    seq_initialized = true;
  }

  const auto& coords = index.at(id);
  const auto seq_len = coords.seq_len;
  btllib::check_error(seq_len >= max_seqlen,
                      FN_NAME + ": Seq size over max limit.");
  seqs_file->seekg(coords.seq_start);
  seqs_file->read(seq, seq_len);
  seq[seq_len] = '\0';

  return { seq, seq_len };
}

#endif