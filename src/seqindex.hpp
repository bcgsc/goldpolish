#ifndef SEQINDEX_HPP
#define SEQINDEX_HPP

#include "utils.hpp"

#include "btllib/status.hpp"

#include <string>
#include <unordered_map>
#include <utility>
#include <fstream>

struct SeqCoordinates
{
  size_t seq_start, seq_len;

  SeqCoordinates(const size_t seq_start, const size_t seq_len)
    : seq_start(seq_start)
    , seq_len(seq_len)
  {}
};

class SeqIndex
{

public:
  SeqIndex(const std::string& seqs_filepath);
  SeqIndex(const std::string& index_filepath, const std::string& seqs_filepath);

  SeqIndex(const SeqIndex&) = delete;
  SeqIndex& operator=(const SeqIndex&) = delete;

  void save(const std::string& filepath);

  template<int i>
  std::tuple<const char*, size_t> get_seq(const std::string& id) const;

  size_t get_seq_len(const std::string& id) const;

  bool seq_exists(const std::string& id) const;

private:
  std::string seqs_filepath;
  std::unordered_map<std::string, SeqCoordinates> seqs_coords;
};

template<int i>
std::tuple<const char*, size_t>
SeqIndex::get_seq(const std::string& id) const
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

  const auto& coords = seqs_coords.at(id);
  const auto seq_len = coords.seq_len;
  btllib::check_error(seq_len >= max_seqlen,
                      FN_NAME + ": Seq size over max limit.");
  seqs_file->seekg(coords.seq_start);
  seqs_file->read(seq, seq_len);
  seq[seq_len] = '\0';

  return { seq, seq_len };
}

#endif