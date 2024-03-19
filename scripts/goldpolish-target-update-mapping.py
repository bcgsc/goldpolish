"""Updates verbose mapping file based on gap sequences and their flanks """
# imports
import argparse
import csv
from intervaltree import Interval, IntervalTree
import btllib


class PafFileRow:
    """Represents row found in paf file"""
    def __init__(self, paf_file):
        self.query_name = paf_file[0]
        self.query_len = int(paf_file[1])
        self.query_start = int(paf_file[2])
        self.query_end = int(paf_file[3])
        self.relative_strand = paf_file[4]
        self.target_name = paf_file[5]
        self.target_length = int(paf_file[6])
        self.target_start = int(paf_file[7])
        self.target_end = int(paf_file[8])
        self.num_matches = int(paf_file[9])
        self.aln_length = int(paf_file[10])
        self.quality = int(paf_file[11])

    def row_to_list(self):
        """Returns class object as a list"""
        return [
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.relative_strand,
            self.target_name,
            self.target_length,
            self.target_start,
            self.target_end,
            self.num_matches,
            self.aln_length,
            self.quality,
        ]


def parse_args():
    """Reads input arguments from command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--gaps",
        help="extracted gaps file in fasta format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mapping",
        help="mapping file, either verbose_mapping or paf",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="name of updated mapping file",
        type=str,
        required=True,
        default="GoldPolish-Target_mapping_updated",
    )
    # parser.add_argument("-v", "--version", action='version', version='version1')
    return parser.parse_args()


def parse_description(description):
    """parses description of fasta file and returns corresponding coordinate interval"""
    start = int(description.split("-")[0])
    end = int(description.split("-")[1]) + 1
    return Interval(start, end)


def populate_trees(trees, name, interval):
    """adds intervals to the interval trees in tree_dict"""
    name_index = name.split(".")
    trees[name_index[0]].addi(interval.begin, interval.end, name)

def make_interval_tree(tree_dict, args):
    "Parses fasta sequences"
    with btllib.SeqReader(args.gaps, btllib.SeqReaderFlag.LONG_MODE) as reader:
        for record in reader:
            seq_name, seq_description = record.id, record.comment
            generic_name = seq_name.split(".")[0]
            if generic_name not in tree_dict:
                tree_dict[generic_name] = IntervalTree()
            seq_interval = parse_description(seq_description)
            populate_trees(tree_dict, seq_name, seq_interval)
        return tree_dict

def update_paf_file(args, tree_dict):
    "Updates paf file based on new coordinate system"
    with open(args.mapping, encoding="utf-8") as f_in, open(
        args.output, "w", encoding="utf-8"
    ) as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        reader = csv.reader(f_in, delimiter="\t", quotechar='"')
        # iterates through all reads in the paffile
        for row in reader:
            row_info = PafFileRow(row)

            contig_name = row_info.target_name
            start_pos = int(row_info.target_start)
            end_pos = int(row_info.target_end)
            query_start = int(row_info.query_start)
            query_end = int(row_info.query_end)

            update_row = False

            # if the contig is in the tree dictionary, checks for mapped minimizers within gap
            if contig_name in tree_dict:
                tree = tree_dict[contig_name]
                new_row = PafFileRow(row)
                # if start/end position maps to gap sequence, updates position
                tree_overlap = tree.overlap(start_pos, end_pos)
                if tree_overlap:
                    update_row = True  # update row in df
                    seq_interval = list(tree_overlap)[0]
                    gap_name = list(tree_overlap)[0].data
                    gap_start = seq_interval.begin
                    gap_end = seq_interval.end

                    new_row.target_name = gap_name
                    # query start and end unchanged
                    if gap_start < start_pos and gap_end > end_pos:
                        new_row.target_start = start_pos - gap_start
                        new_row.target_end = end_pos - gap_start
                    # query end updated
                    elif gap_start < start_pos and gap_end <= end_pos:
                        new_row.target_start = start_pos - gap_start
                        new_row.target_end = gap_end - gap_start

                        new_row.query_end = query_end - (end_pos - gap_end)
                    # query start updated
                    elif gap_start >= start_pos and gap_end > end_pos:
                        new_row.target_start = 0
                        new_row.target_end = end_pos - gap_start

                        new_row.query_start = query_start + (gap_start - start_pos)
                    # query start and end updated
                    else:
                        new_row.target_start = 0
                        new_row.target_end = gap_end - gap_start

                        new_row.query_start = query_start + (gap_start - start_pos)
                        new_row.query_end = query_end - (end_pos - gap_end)
                    
                    new_row.aln_length = new_row.target_end - new_row.target_start
                    new_row_list = new_row.row_to_list()

                if update_row:
                    writer.writerow(new_row_list)


def main():
    "Run update mapping file"
    args = parse_args()
    # makes tree dictionary with contig names
    tree_dict = {}
    trees = make_interval_tree(tree_dict, args)

    if "paf" in args.mapping:
        update_paf_file(args, trees)
    else:
        raise ValueError("Mapping file must be a paf file")

if __name__ == "__main__":
    main()
