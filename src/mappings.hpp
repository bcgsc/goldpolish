#ifndef MAPPINGS_HPP
#define MAPPINGS_HPP

#include "seqindex.hpp"

#include <string>
#include <unordered_map>

struct Mapping
{
  std::string seq_id;
  unsigned mx_in_common;

  Mapping(std::string seq_id, const unsigned mx_in_common)
    : seq_id(std::move(seq_id))
    , mx_in_common(mx_in_common)
  {}
};

using Mappings = std::vector<Mapping>;

class AllMappings
{

public:
  AllMappings(const std::string& filepath,
              const SeqIndex& target_seqs_index,
              unsigned mx_threshold_min,
              unsigned mx_threshold_max,
              double mx_max_mapped_seqs_per_target_10kbp);

  AllMappings(const AllMappings&) = delete;
  AllMappings& operator=(const AllMappings&) = delete;

  const Mappings& get_mappings(const std::string& id) const;

private:
  void load_ntlink(const std::string& filepath,
                   const SeqIndex& target_seqs_index,
                   unsigned mx_threshold_min);
  void load_sam(const std::string& filepath, const SeqIndex& target_seqs_index);

  void filter(double max_mapped_seqs_per_target_10kbp,
              unsigned mx_threshold_min,
              unsigned mx_threshold_max,
              const SeqIndex& target_seqs_index);

  std::unordered_map<std::string, Mappings> all_mappings;
  static const Mappings EMPTY_MAPPINGS;
};

#endif
