"""Extracts regions to polish with additional flanks into fasta file"""
# imports
import argparse
import csv
import re
import btllib
from collections import namedtuple

Coordinate = namedtuple("Coordinate", "start end")
min_gap_length = 1

def parse_args():
    """Parses arguments passed in command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", help="target file in fasta format", type=str, required=True
    )
    parser.add_argument(
        "-b",
        "--bed",
        help="bed file specifying regions to polish",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="prefix of output file [<output>.fa]",
        type=str,
        required=False,
        default="GoldPolish-Target_extracted_gaps.fa",
    )
    parser.add_argument(
        "-l",
        "--length",
        help="length of flanking regions [64]",
        type=int,
        required=False,
        default=64,
    )
    return parser.parse_args()


def make_coord_dict(bed):
    """Creates a dictionary of tuples representing regions to polish"""
    coord_dict_2 = {}

    with open(bed, encoding="utf-8") as bed_file:
        bed_reader = csv.reader(bed_file, delimiter="\t", quotechar='"')
        for row in bed_reader:
            contig_name = row[0]
            coord = Coordinate(row[1], row[2])

            if contig_name not in coord_dict_2:
                coord_dict_2[contig_name] = [coord]
            else:
                coord_dict_2[contig_name].append(coord)
    return coord_dict_2


def extract_masked_subsequences(sequence, name, flank_length, writer):
    """extracts lowercase sequences and returns list of sequences with flank_length bp long flanking regions"""
    gap_count = 1
    idx = 1  # index of subseqs

    subseqs = re.findall(r"([A-Z]+|[a-z]+)", sequence)
    filtered_subseqs = [] # short uppercase seqs are appended to adjacent softmasked seqs
    filtered_subseqs.append(
        subseqs[0]
    )  # first item always appended, avoids 0-index issues

    while idx < len(subseqs):
        subseq = subseqs[idx]
        if subseq.isupper():
            if len(subseq) >= 2 * flank_length:  # surpases min length threshold
                filtered_subseqs.append(subseq)
            elif filtered_subseqs[-1].islower():
                # append lower case subseq to previous subseq
                filtered_subseqs[-1] = filtered_subseqs[-1] + subseq.lower()
            else:
                raise Exception("Unexpected order of softmasked subseqs")
        else:  # subseq is lower
            if filtered_subseqs[-1].isupper():
                filtered_subseqs.append(subseq)
            elif filtered_subseqs[-1].islower():
                # previous subseq is lower, need to concat with prev
                filtered_subseqs[-1] = filtered_subseqs[-1] + subseq
            else:
                raise Exception("Unexpected order of softmasked subseqs")
        idx += 1
        
    idx = 0

    for subseq in filtered_subseqs:
        if subseq.islower() and len(subseq) > min_gap_length:
            flank_start = max(0, idx - flank_length)
            flank_end = min(len(sequence), idx + len(subseq) + flank_length)
            if flank_end >= flank_start:  # just in case
                flanked_subseq = sequence[flank_start:flank_end]

                write_flanked_subsequence(
                    flanked_subseq,
                    flank_start,
                    flank_end,
                    name + "." + str(gap_count),
                    writer,
                )
                gap_count += 1
        idx += len(subseq)


def extract_subsequences_from_bed(sequence, name, flank_length, writer, coords):
    """extracts seqs with coordinates from bed file and returns sequences with flank_length bp long flanks"""
    if name in coords:
        count = 0
        coord_list = coords[name]
        idx = 1

        filtered_coords = [] # short uppercase seqs appended to adj coordinates
        filtered_coords.append(
            coord_list[0]
        )  # first item always appended, avoids 0-index issues

        while idx < len(coord_list):
            coord = Coordinate(*coord_list[idx])
            prev_coord = Coordinate(*filtered_coords[-1])

            if (
                int(coord.start) - int(prev_coord.end)
            ) < 2 * flank_length:  # length between adjacent coords too small
                filtered_coords[-1] = (prev_coord.start, coord.end)
            else:
                filtered_coords.append(coord_list[idx])
            idx += 1

        for coord in filtered_coords:
            start = max(0, int(coord[0]) - flank_length)
            end = min(len(sequence), int(coord[1]) + flank_length)

            count += 1

            write_flanked_subsequence(
                sequence[start:end], start, end, name + "." + str(count), writer
            )


def write_flanked_subsequence(subsequence, start_flank, end_flank, gap_name, writer):
    """Writes subsequence with flanking regions into fasta file"""
    writer.write(gap_name, str(start_flank) + "-" + str(end_flank), subsequence.upper())


def main():
    """Parses fasta file to extract sequences + flanks"""
    args = parse_args()
    writer_fasta = btllib.SeqWriter(args.output, btllib.SeqWriter.FASTA)

    # makes coordinate dictionary if bed file provided
    if args.bed:
        coord_dict = make_coord_dict(args.bed)

    # loop through sequences in fasta file
    with btllib.SeqReader(args.fasta, btllib.SeqReaderFlag.LONG_MODE) as reader:
        for record in reader:
            seq_name, seq = record.id, str(record.seq)
            if not args.bed:
                extract_masked_subsequences(
                    seq, seq_name, int(args.length), writer_fasta
                )
            else:
                if seq_name in coord_dict:
                    extract_subsequences_from_bed(
                        seq, seq_name, int(args.length), writer_fasta, coord_dict
                    )

    writer_fasta.close()

if __name__ == "__main__":
    main()
