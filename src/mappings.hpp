#ifndef MAPPINGS_HPP
#define MAPPINGS_HPP

#include "seqindex.hpp"

#include <set>
#include <string>
#include <unordered_map>

using SeqId = std::string;

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

  const std::vector<SeqId>& get_mappings(const std::string& id) const;

private:
  void load_ntlink(const std::string& filepath,
                   const SeqIndex& target_seqs_index,
                   unsigned mx_threshold_min);
  void load_sam(const std::string& filepath, const SeqIndex& target_seqs_index);
  void load_paf(const std::string& filepath, const SeqIndex& target_seqs_index);

  void filter(double max_mapped_seqs_per_target_10kbp,
              unsigned mx_threshold_min,
              unsigned mx_threshold_max,
              const SeqIndex& target_seqs_index);

  void load_mapping(const std::string& mapped_seq_id,
                    const std::string& target_seq_id,
                    const SeqIndex& target_seqs_index,
                    unsigned mx = 0);

  std::unordered_map<SeqId, std::vector<SeqId>> all_mappings;
  std::unordered_map<SeqId, std::set<SeqId>> all_inserted_mappings;
  std::unordered_map<SeqId, std::vector<unsigned>> all_mx_in_common;

  static const std::vector<SeqId> EMPTY_MAPPINGS;
};

#endif
