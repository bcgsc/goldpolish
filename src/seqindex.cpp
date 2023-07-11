#include "seqindex.hpp"
#include "utils.hpp"

#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

SeqIndex::SeqIndex(const std::string& seqs_filepath)
  : seqs_filepath(seqs_filepath)
{
  btllib::log_info(FN_NAME + ": Building index for " + seqs_filepath + "... ");

  std::ifstream seqsfile(seqs_filepath);
  btllib::check_stream(seqsfile, seqs_filepath);

  const auto fastq = (seqsfile.peek() == '@');

  std::string line;
  std::string id;
  long i = 0, byte = 0 /*, id_startbyte = 0*/, id_endbyte = 0, seq_start = 0,
       seq_len = 0;
  while (bool(std::getline(seqsfile, line))) {
    const long endbyte = byte + long(line.size());
    if (fastq) {
      if (i % 4 == 0) {
        /*id_startbyte = byte + 1;*/
        id_endbyte = endbyte;
        id = btllib::split(line, " ")[0].substr(1);
      } else if (i % 4 == 1) {
        seq_start = id_endbyte + 1;
        seq_len = endbyte - id_endbyte - 1;

      } else if (i % 4 == 3) {

        seqs_coords_and_phred_avg.emplace(
          std::piecewise_construct,
          std::make_tuple(id),
          std::make_tuple(seq_start,
                          seq_len,
                          btllib::calc_phred_avg(line, 0, line.size() - 1)));
      }
    } else {
      if (i % 2 == 0) {
        /*id_startbyte = byte + 1;*/
        id_endbyte = endbyte;
        id = btllib::split(line, " ")[0].substr(1);
      } else {
        const auto seq_start = id_endbyte + 1;
        const auto seq_len = endbyte - id_endbyte - 1;
        seqs_coords_and_phred_avg.emplace(
          std::piecewise_construct,
          std::make_tuple(id),
          std::make_tuple(seq_start, seq_len, 0.0));
      }
    }
    byte = endbyte + 1;
    i++;
  }

  btllib::log_info(FN_NAME + ": Done.");
}

void
SeqIndex::save(const std::string& filepath)
{
  btllib::log_info(FN_NAME + ": Saving index to " + filepath + "... ");

  std::ofstream indexfile(filepath);
  for (const auto& seq_coords : seqs_coords_and_phred_avg) {
    const auto id = seq_coords.first;
    const auto seq_start = seq_coords.second.seq_start;
    const auto seq_len = seq_coords.second.seq_len;
    const auto phred_avg = seq_coords.second.phred_avg;
    indexfile << id << '\t' << seq_start << '\t' << seq_len << '\t' << phred_avg
              << '\n';
  }

  btllib::log_info(FN_NAME + ": Done.");
}

SeqIndex::SeqIndex(const std::string& index_filepath, std::string seqs_filepath)
  : seqs_filepath(std::move(seqs_filepath))
{
  btllib::log_info(FN_NAME + ": Loading index from " + index_filepath + "... ");

  std::ifstream ifs(index_filepath);
  btllib::check_stream(ifs, index_filepath);
  std::string token, id;
  unsigned long seq_start = 0, seq_len = 0;
  double phred_avg = 0.0;
  unsigned long i = 0;
  while (bool(ifs >> token)) {
    switch (i % 4) {
      case 0:
        id = std::move(token);
        break;
      case 1:
        seq_start = std::stoul(token);
        break;
      case 2: {
        seq_len = std::stoul(token);
        break;
      }
      case 3: {
        phred_avg = std::stod(token);
        seqs_coords_and_phred_avg.emplace(
          std::piecewise_construct,
          std::make_tuple(id),
          std::make_tuple(seq_start, seq_len, phred_avg));
        break;
      }
      default: {
        btllib::log_error(FN_NAME + ": Invalid switch branch.");
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }
    i++;
  }
  btllib::log_info(FN_NAME + ": Done!");
}

size_t
SeqIndex::get_seq_len(const std::string& id) const
{
  return seqs_coords_and_phred_avg.at(id).seq_len;
}

double
SeqIndex::get_phred_avg(const std::string& id) const
{
  return seqs_coords_and_phred_avg.at(id).phred_avg;
}

bool
SeqIndex::seq_exists(const std::string& id) const
{
  return seqs_coords_and_phred_avg.find(id) != seqs_coords_and_phred_avg.end();
}