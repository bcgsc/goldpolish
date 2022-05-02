#ifndef UTILS_HPP
#define UTILS_HPP

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

using SeqId = std::string;

struct ReadMapping {
  SeqId seq_id;
  unsigned mx_in_common;

  ReadMapping(const SeqId& seq_id, const unsigned mx_in_common) : seq_id(seq_id), mx_in_common(mx_in_common) {}
};

struct SeqCoordinates {
  size_t seq_start, seq_len;

  SeqCoordinates(const size_t seq_start, const size_t seq_len)
    : seq_start(seq_start)
    , seq_len(seq_len) {}
};

using ReadsMapping = std::unordered_map<SeqId, std::vector<ReadMapping>>;
using Index = std::unordered_map<SeqId, SeqCoordinates>;

std::vector<size_t>
get_random_indices(size_t total_size, size_t count);

void
load_index(Index& index, const std::string& filepath);

void
load_reads_mapping(ReadsMapping& reads_mapping, const std::string& filepath, const Index& contigs_index, unsigned mx_threshold_min);

void
filter_read_mappings(ReadsMapping& reads_mapping, double max_reads_per_contig_10kbp, unsigned mx_threshold_min, unsigned mx_threshold_max, Index& contigs_index);

void
load_reads_mapping_sam(ReadsMapping& reads_mapping, const std::string& filepath, const Index& contigs_index);

void
fill_bfs(const char* seq,
         size_t seq_len,
         unsigned hash_num,
         const std::vector<unsigned>& ks,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs);

inline void
fill_bfs(const std::string& seq,
         unsigned hash_num,
         const std::vector<unsigned>& ks,
         unsigned kmer_threshold,
         std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>>& cbfs,
         std::vector<std::unique_ptr<btllib::KmerBloomFilter>>& bfs) {
   fill_bfs(seq.c_str(), seq.size(), hash_num, ks, kmer_threshold, cbfs, bfs);
}

void
bind_to_parent();

template<int i>
std::tuple<const char*, size_t>
get_seq_with_index(const std::string& id,
                   const Index& index,
                   const std::string& seqs_filepath)
{
  thread_local static std::ifstream *seqs_file;
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
  btllib::check_error(seq_len >= max_seqlen, "Read size over max.");
  seqs_file->seekg(coords.seq_start);
  seqs_file->read(seq, seq_len);
  seq[seq_len] = '\0';

  return { seq, seq_len };
}

#endif