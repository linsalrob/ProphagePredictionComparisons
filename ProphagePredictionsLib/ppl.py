"""
Some common methods used in this application
"""

import binascii
import gzip
from Bio import SeqIO


def is_gzip(gbkf):
    """
    Is the file compressed?
    :param gbkf:
    :return: true if compressed else false
    """

    with open(gbkf, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def genbank_seqio(gbkf):
    """
    Get the parser stream
    :param gbkf: genbank file
    :return:
    """

    if is_gzip(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')
    return SeqIO.parse(handle, "genbank")


def feature_id(seq, feat):
    """
    Choose the appropriate id for the feature
    :param seq: the bioseq object
    :param feat: the feature
    :return: the id
    """

    if 'protein_id' in feat.qualifiers:
        return '|'.join(feat.qualifiers['protein_id'])
    elif 'locus_tag' in feat.qualifiers:
        return "|".join(feat.qualifiers['locus_tag'])
    elif 'db_xref' in feat.qualifiers:
        return '|'.join(feat.qualifiers['db_xref'])
    else:
        return seq.id + "." + str(feat.location)
