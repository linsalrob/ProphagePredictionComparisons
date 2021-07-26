from Bio import SeqIO
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-in', '--input', help='genbank file', required=True)
    parser.add_argument('-out', '--output', help='genbank file with simplified IDs')
    args = parser.parse_args()

    seqs = SeqIO.parse(args.input, "genbank")
    sanitized_seqs = []
    for seq in seqs:
        if '.' in seq.id:
            seq.id = seq.name = seq.id.replace('.', '_dot_')
        sanitized_seqs.append(seq)
    SeqIO.write(sanitized_seqs, args.output, "genbank")