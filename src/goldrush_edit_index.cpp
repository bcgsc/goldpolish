#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

void
build_index(const std::string& seqs_filepath, const std::string& index_filepath)
{
  std::ifstream seqsfile(seqs_filepath);
  std::ofstream indexfile(index_filepath);

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
        indexfile << id << '\t' << id_startbyte << '\t' << id_endbyte << '\t'
                  << endbyte << '\n';
      }
    } else {
      if (i % 2 == 0) {
        id_startbyte = byte + 1;
        id_endbyte = endbyte;
        id = btllib::split(line, " ")[0].substr(1);
      } else {
        indexfile << id << '\t' << id_startbyte << '\t' << id_endbyte << '\t'
                  << endbyte << '\n';
      }
    }
    byte = endbyte + 1;
    i++;
  }
}

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

  build_index(seqs_filepath, index_filepath);

  return 0;
}