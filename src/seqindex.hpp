#ifndef SEQINDEX_HPP
#define SEQINDEX_HPP

#include "utils.hpp"

#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>

#include <fcntl.h>
#include <unistd.h>

struct SeqCoordinatesAndPhredAvg
{
  size_t seq_start, seq_len;
  double phred_avg;

  SeqCoordinatesAndPhredAvg(const size_t seq_start,
                            const size_t seq_len,
                            const double phred_avg)
    : seq_start(seq_start)
    , seq_len(seq_len)
    , phred_avg(phred_avg)
  {
  }
};

class SeqIndex
{

public:
  SeqIndex(const std::string& seqs_filepath);
  SeqIndex(const std::string& index_filepath, std::string seqs_filepath);

  SeqIndex(const SeqIndex&) = delete;
  SeqIndex& operator=(const SeqIndex&) = delete;

  void save(const std::string& filepath);

  template<int i>
  std::tuple<const char*, size_t> get_seq(const std::string& id) const;

  size_t get_seq_len(const std::string& id) const;

  double get_phred_avg(const std::string& id) const;

  bool seq_exists(const std::string& id) const;

private:
  std::string seqs_filepath;
  std::unordered_map<std::string, SeqCoordinatesAndPhredAvg>
    seqs_coords_and_phred_avg;
};

template<int i>
std::tuple<const char*, size_t>
SeqIndex::get_seq(const std::string& id) const
{
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  thread_local static int seqs_file;
  thread_local static bool seqs_file_initialized = false;

  static const size_t max_seqlen = 20ULL * 1024ULL * 1024ULL;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  thread_local static char* seq;
  thread_local static bool seq_initialized = false;

  if (!seqs_file_initialized) {
    seqs_file = open(seqs_filepath.c_str(), O_RDONLY);
#ifdef __linux__
    const auto ret = posix_fadvise(seqs_file, 0, 0, POSIX_FADV_RANDOM);
    btllib::check_error(ret != 0, "posix_fadvise failed.");
#endif
    seqs_file_initialized = true;
  }

  if (!seq_initialized) {
    seq = new char[max_seqlen];
    seq_initialized = true;
  }

  const auto& coords = seqs_coords_and_phred_avg.at(id);
  const auto seq_len = coords.seq_len;
  btllib::check_error(
    seq_len >= max_seqlen,
    FN_NAME + ": Seq size over buffer size. Consider increasing buffer size.");
  // NOLINTNEXTLINE(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
  const auto lseek_ret = lseek(seqs_file, coords.seq_start, SEEK_SET);
  btllib::check_error(lseek_ret == -1,
                      FN_NAME + ": lseek: " + btllib::get_strerror());
  const auto read_ret = read(seqs_file, seq, seq_len);
  btllib::check_error(read_ret != int(seq_len),
                      FN_NAME + ": read did not read all bytes.");

  seq[seq_len] = '\0';

  return decltype(SeqIndex::get_seq<i>(id)){ seq, seq_len };
}

#endif
