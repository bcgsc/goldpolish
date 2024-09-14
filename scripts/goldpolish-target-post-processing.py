"""Inserts polished gap sequences back into the original assembly"""
# imports
import argparse
import btllib


def parse_args():
    """parses args passed through command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fasta", help="target file in fasta format", type=str, required=True
    )
    parser.add_argument(
        "-g",
        "--gaps",
        help="polished gaps and flanks in fasta format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="name of output file",
        type=str,
        required=True
    )
    return parser.parse_args()


def insert_seq(gap_coords, sequence):
    """inserts gap_sequences into contig at the given coordinates"""
    start = 0
    updated_seq = []

    # reading gap names to get coordinates
    for gap in gap_coords:
        gap_start = gap[1]
        gap_end = gap[2]
        gap_sequence = gap[3]

        updated_seq.append(sequence[start : int(gap_start)])
        updated_seq.append(gap_sequence)
        start = int(gap_end)

    # just to check
    if start < len(sequence) - 1:
        updated_seq.append(sequence[int(start):])
    return "".join(updated_seq)


def make_seq_dict(gaps):
    with btllib.SeqReader(gaps, btllib.SeqReaderFlag.LONG_MODE) as gap_reader:
        gap_count_dict = {}
        for record in gap_reader:
            gap_name = record.id
            gap_name_info = gap_name.rsplit(".", 1)
            gap_count_dict[".".join(gap_name_info[:-1])] = gap_name_info[-1]
    return gap_count_dict


def insert_and_write_seqs(fasta, gaps, gap_dict, writer):
    "Reinserts sequences into the fasta file using gap coordinates"
    with btllib.SeqReader(fasta, btllib.SeqReaderFlag.LONG_MODE) as reader:
        with btllib.SeqReader(gaps, btllib.SeqReaderFlag.LONG_MODE) as gap_reader:
            gap_record = gap_reader.read()
            for record in reader:
                name, sequence = record.id, str(record.seq)
                if name not in gap_dict:
                    writer.write(name, "", sequence)
                    continue
                gap_name = gap_record.id.rsplit(".", 1)[0] + ".0"
                gap_coordinates = []

                # making sure you don't increment past end of sequence with gap_reader.read()
                while gap_name.rsplit(".", 1)[0] == name and int(
                    gap_name.rsplit(".", 1)[-1]
                ) < int(gap_dict[name]):
                    gap_name = gap_record.id
                    gap_desc = gap_record.comment
                    gap_seq = gap_record.seq
                    start = gap_desc.split("-")[0]
                    end = gap_desc.split("-")[1]

                    gap_coordinates.append(
                        (
                            gap_name,
                            start,
                            end,
                            gap_seq,
                        )
                    )
                    gap_record = gap_reader.read()

                sequence = insert_seq(gap_coordinates, sequence)
                writer.write(name, "", sequence)


def main():
    # setting up inputs and outputs
    args = parse_args()

    writer = btllib.SeqWriter(args.output, btllib.SeqWriter.FASTA)
    gap_count_dict = make_seq_dict(args.gaps)

    insert_and_write_seqs(args.fasta, args.gaps, gap_count_dict, writer)
    writer.close()

if __name__ == "__main__":
    main()