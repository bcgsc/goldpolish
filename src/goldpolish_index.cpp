#include "seqindex.hpp"

int
main(int argc, char** argv)
{
  if (argc != 3) {
    std::cerr << "Wrong args.\n";
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  unsigned arg = 1;
  auto* const seqs_filepath = argv[arg++];
  auto* const index_filepath = argv[arg++];

  SeqIndex index(seqs_filepath);
  index.save(index_filepath);

  return 0;
}