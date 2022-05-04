#include "seqindex.hpp"

int
main(int argc, char** argv)
{
  if (argc != 3) {
    std::cerr << "Wrong args.\n";
    std::exit(EXIT_FAILURE);
  }
  unsigned arg = 1;
  const auto seqs_filepath = argv[arg++];
  const auto index_filepath = argv[arg++];

  SeqIndex index(seqs_filepath);
  index.save(index_filepath);

  return 0;
}