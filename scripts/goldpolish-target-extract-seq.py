"""Extracts regions to polish with additional flanks into fasta file"""
# imports
import argparse
import csv
import re
import btllib


def parse_args():
    """Parses arguments passed in command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", help="target file in fasta format", type=str, required=True
    )
    parser.add_argument(
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
    # parser.add_argument("-v", "--version", action='version', version='version1')
    return parser.parse_args()


def make_coord_dict(bed):
    """Creates a dictionary of tuples representing regions to polish"""
    coord_dict_2 = {}

    with open(bed, encoding="utf-8") as bed_file:
        bed_reader = csv.reader(bed_file, delimiter="\t", quotechar='"')
        for row in bed_reader:
            contig_name = row[0]
            start = row[1]
            end = row[2]

            if contig_name not in coord_dict_2:
                coords = []

            coords.append((start, end))
            coord_dict_2[contig_name] = coords
    return coord_dict_2


def extract_masked_subsequences(sequence, name, length, writer):
    """extracts lowercase sequences and returns list of sequences with l bp long flanking regions"""
    gap_count = 1
    idx = 1

    subseqs = re.findall(r"([A-Z]+|[a-z]+)", sequence)

    while idx < len(subseqs):
        subseq = subseqs[idx]
        if subseq and subseq.isupper() and len(subseq) < 2 * length:
            first_subseq = subseqs[idx - 1] if idx > 0 and subseqs[idx - 1] else ""
            second_subseq = (
                subseqs[idx + 1] if idx < len(subseqs) - 1 and subseqs[idx + 1] else ""
            )
            updated_subseq = first_subseq + subseq.lower() + second_subseq

            subseqs = (
                subseqs[: idx - 1]
                + [None]
                + [None]
                + [updated_subseq]
                + subseqs[idx + 2 :]
            )
            idx += 1
        idx += 1
    # getting rid of old indices added in last step
    subseqs = filter(None, subseqs)
    idx = 0

    for subseq in subseqs:
        if subseq.islower() and len(subseq) > 1:
            flank_start = max(0, idx - length)
            flank_end = min(len(sequence), idx + len(subseq) + length)
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


def extract_subsequences_from_bed(sequence, name, length, writer, coords):
    """extracts seqs with coordinates from bed file and returns sequences with l bp long flanks"""
    if name in coords:
        count = 0
        coord_list = coords[name]
        idx = 0

        while idx < len(coord_list) - 1:
            # distance between coords of seqs to be extracted
            inter_seq_len = int(coord_list[idx + 1][0]) - int(coord_list[idx][1])
            if inter_seq_len < 2 * length:
                start = coord_list[idx][0]
                end = coord_list[idx + 1][1]

                list_before_index = coord_list[:idx] if idx > 0 else []
                list_after_index = (
                    coord_list[idx + 2 :] if idx < len(coord_list) + 2 else []
                )
                coord_list = (
                    list_before_index + [None] + [(start, end)] + list_after_index
                )
            idx += 1

        coord_list = filter(None, coord_list)

        for coord in coord_list:
            start = max(0, int(coord[0]) - length)
            end = min(len(sequence), int(coord[1]) + length)

            count += 1

            write_flanked_subsequence(
                sequence[start:end], start, end, name + "." + str(count), writer
            )


def write_flanked_subsequence(subsequence, start_flank, end_flank, gap_name, writer):
    """Writes subsequence with flanking regions into fasta file"""
    writer.write(gap_name, str(start_flank) + "-" + str(end_flank), subsequence.upper())


def main():
    "Parses fasta file to extract sequences + flanks"
    args = parse_args()
    writer_fasta = btllib.SeqWriter(args.output, btllib.SeqWriter.FASTA)

    # makes coordinate dictionary if bed file provided
    if args.bed != "":
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