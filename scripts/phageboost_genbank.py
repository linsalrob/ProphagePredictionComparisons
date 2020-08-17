"""
Run phage boost.

This is extracted from their example code, but we added parsing genbank files to dataframes

"""

import sys
import argparse

import pandas as pd
from PhageBoost.main import calculate_features, read_model_from_file, predict, get_predictions

from ProphagePredictionsLib import genbank_seqio, feature_id

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def genbank_to_pandas(gbkf, mincontiglen, ignorepartials=True, convert_selenocysteine=False):
    """
    This is a bit of a specific format used by phage_boost. its a simple dataframe with a couple of
    additional columns:
        ['contig',
         'id',
         'start',
         'stop',
         'direction',
         'partial',
         'DNAseq',
         'AAseq',
         'header']
    :param mincontiglen: minimum contig length
    :param ignorepartials: Ignore any gene call with a frameshift (ie. a stop codon in the middle of the sequence)
    :param convert_selenocysteine: PhageBoost crashes with a selenocysteine protein because it is not in Biopython
    :param gbkf: Genbank file to parse
    :return: a pandas data frame
    """

    c = 0
    genes = []
    for seq in genbank_seqio(gbkf):
        if len(seq) < mincontiglen:
            sys.stderr.write(f"Skipped {seq.id} because it's length ({len(seq)}) is less than the " +
                             "minimum contig length ({mincontiglen})\n")
            continue
        for feat in seq.features:
            if feat.type != 'CDS':
                continue

            tid = seq.id + "_" + str(c)
            partial = 0
            # I don't think this is exactly right
            if 'truncated' in feat.qualifiers:
                partial = 1

            dnaseq = str(feat.extract(seq).seq)
            if len(dnaseq) == 0:
                sys.stderr.write(f"The DNA sequence for {feature_id(seq, feat)} was zero, so skipped\n")
                continue

            # we just do a de novo translation rather than relying on the translation provided
            # in the genbank file that is often wrong
            trans = str(feat.extract(seq).translate().seq)

            while trans.endswith('*'):
                trans = trans[:-1]

            # Partial amino acid codes we should ignore. These are not present in BioPython's SeqUtils::ProtParam
            # and it crashes the system
            paa = {'B', 'Z', 'J', 'X', '*'}

            keeporf = True

            if ignorepartials:
                for aa in paa:
                    if aa in trans:
                        sys.stderr.write(f"There is a {aa} in  {feature_id(seq, feat)} so skipped.\n")
                        keeporf = False

            if not keeporf:
                continue

            if len(trans) == 0:
                sys.stderr.write(f"The translation for {feature_id(seq, feat)} was zero, so skipped.\n")
                continue

            if convert_selenocysteine:
                trans = trans.replace('U', 'C')
            row = [seq.id, c, feat.location.start.position, feat.location.end.position, feat.strand,
                   partial, dnaseq, trans, tid]
            c += 1

            genes.append(row)

    gc = pd.DataFrame(genes, columns=['contig', 'id', 'start', 'stop', 'direction', 'partial', 'DNAseq', 'AAseq',
                                      'header'])
    return gc


def run_phage_boost(gcs, model_file, verbose):
    """
    Run phage boost
    :param model_file: The model file that is probably something like model_delta_std_hacked.pickled.silent.gz
    :param gcs: The pandas data frame of gene calls
    :param verbose: more output
    :return:
    """
    # rolling params
    period = 20
    win_type = 'parzen'
    min_periods = 1

    # region finding params
    threshold = 0.9
    length = 10
    gaps = 5
    neighbouring = 0
    alpha = 0.001

    # calculate features from gene calls
    if verbose:
        sys.stderr.write("Calculating features\n")

    df = calculate_features(gcs)
    # load model
    model, feats, feats_, limit = read_model_from_file(model_file)
    # transform data
    df = get_predictions.get_deltas(df[feats_])
    if verbose:
        sys.stderr.write("Transforming gene predictions to regions\n")
    # transform single gene predictions to regions
    newgenecalls, nphages, r = predict(model, gcs, df, feats, period, win_type,
                                       min_periods, limit, threshold, length,
                                       gaps, neighbouring, alpha)
    return r


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', '--genbankfile', help='GenBank file to parse', required=True)
    parser.add_argument('-m', '--modelfile', required=True,
                        help="Model file. Probably something like model_delta_std_hacked.pickled.silent.gz")
    parser.add_argument('-o', '--outputfile', help='output file for phage regions')
    parser.add_argument('-c', '--mincontiglen', default=1000, type=int,
                        help='minimum contig length  [Default: %(default)d]')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        sys.stderr.write("Reading genbank file\n")
    genecalls = genbank_to_pandas(args.genbankfile, args.mincontiglen, True, True)
    if args.verbose:
        sys.stderr.write("Phage Boosting\n")
    res = run_phage_boost(genecalls, args.modelfile, args.verbose)
    if args.outputfile:
        with open(args.outputfile, 'w') as out:
            res.to_csv(out, sep="\t", header=True)
    else:
        print(res)
