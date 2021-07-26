from seeker import SeekerFasta
import argparse


def run_seeker(input, output, opt_fasta, lstm_type):
    seeker_fasta = SeekerFasta(f"{input}", LSTM_type=lstm_type)
    seeker_fasta.save2bed(f"{output}")
    if opt_fasta:
        seeker_fasta.save2fasta(f"{output}.fasta")
    print(f"[Seeker] Done {input}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Better way to run Seeker")
    parser.add_argument('-i', '--input', required=True,
                        help="FASTA file to process")
    parser.add_argument('-o', '--output', required=True,
                        help="Directory with output filename (BED file with extension)")
    parser.add_argument('-f', '--fasta', required=False, action='store_false', default=False,
                        help="If provided, FASTA file with found viral sequences will be created")
    parser.add_argument('--lstm', required=False, action='store_const', const="prophage", default="prophage",
                        help="Type of LSTM to be used. 'propage' type is default.")
    args = parser.parse_args()

    run_seeker(args.input, args.output, args.fasta, args.lstm)
