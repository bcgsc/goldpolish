#include "seqindex.hpp"

#include "btllib/status.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

SeqIndex::SeqIndex(const std::string& seqs_filepath)
  : seqs_filepath(seqs_filepath)
{
  btllib::log_info(FN_NAME + ": Building index for " + seqs_filepath + "... ");

  std::ifstream seqsfile(seqs_filepath);

  bool fastq = seqsfile.peek() == '@' ? true : false;

  std::string line;
  std::string id;
  long i = 0, byte = 0, id_startbyte = 0, id_endbyte = 0;
  while (std::getline(seqsfile, line)) {
    const auto endbyte = byte + line.size();
    if (fastq) {
      if (i % 4 == 0) {
        id_startbyte = byte + 1;
        id_endbyte = endbyte;
        id = btllib::split(line, " ")[0].substr(1);
      } else if (i % 4 == 1) {
        const auto seq_start = id_endbyte + 1;
        const auto seq_len = endbyte - id_endbyte - 1;
        seqs_coords.emplace(std::piecewise_construct,
                            std::make_tuple(id),
                            std::make_tuple(seq_start, seq_len));
      }
    } else {
      if (i % 2 == 0) {
        id_startbyte = byte + 1;
        id_endbyte = endbyte;
        id = btllib::split(line, " ")[0].substr(1);
      } else {
        const auto seq_start = id_endbyte + 1;
        const auto seq_len = endbyte - id_endbyte - 1;
        seqs_coords.emplace(std::piecewise_construct,
                            std::make_tuple(id),
                            std::make_tuple(seq_start, seq_len));
      }
    }
    byte = endbyte + 1;
    i++;
  }

  btllib::log_info(FN_NAME + ": Done.");
}

SeqIndex::save(const std::string& filepath)
{
  btllib::log_info(FN_NAME + ": Saving index to " + filepath + "... ");

  std::ofstream indexfile(index_filepath);
  for (const auto& seq_coords : seqs_coords) {
    const auto id = seq_coords.first;
    const auto seq_start = seq_coords.second.seq_start;
    const auto seq_len = seq_coords.second.seq_len;
    indexfile << id << '\t' << seq_start << '\t' << seq_len << '\n';
  }

  btllib::log_info(FN_NAME + ": Done.");
}

SeqIndex::SeqIndex(const std::string& index_filepath,
                   const std::string& seqs_filepath)
  : seqs_filepath(seqs_filepath)
{
  btllib::log_info(FN_NAME + ": Loading index from " + index_filepath + "... ");

  std::ifstream ifs(index_filepath);
  std::string token, id;
  unsigned long seq_start = 0, seq_len = 0;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 4) {
      case 0:
        id = std::move(token);
        break;
      case 1:
        seq_start = std::stoul(token);
        break;
      case 3: {
        seq_len = std::stoul(token);
        seqs_coords.emplace(std::piecewise_construct,
                                  std::make_tuple(id),
                                  std::make_tuple(seq_start, seq_len));
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