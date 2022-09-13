#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Missing pipe path\n";
    return -1;
  }
  std::ofstream of(argv[1]);
  of << 1 << std::flush;
  return 0;
}
