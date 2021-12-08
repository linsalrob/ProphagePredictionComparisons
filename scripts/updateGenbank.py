import sys
import gzip
from Bio import SeqIO

# usage:
# python updateGenbank.py regions.bed genome.gbff.gz > newGenome.gbff
# python updateGenbank.py regions.bed genome.gbff.gz | gzip - > newGenome.gbff.gz

bed = {}
# read in BED coords
with open(sys.argv[1], 'r') as b:
    for line in b:
        if line.startswith('#'):
            continue
        l = line.split()
        if len(l) >=3:
            try:
                bed[l[0]]
            except KeyError:
                bed[l[0]] = []
            bed[l[0]].append([int(l[1]),int(l[2])])

# parse the stupid genbank file
fh = gzip.open(sys.argv[2], 'rt')
genbankfilesarestupid = SeqIO.parse(fh, 'genbank')
for seq in genbankfilesarestupid:
    for i in range(len(seq.features)):
        try:
            del seq.features[i].qualifiers['is_phage']
        except KeyError:
            pass
        try:
            for b in bed[seq.id]:
                if seq.features[i].location.start < b[1] and seq.features[i].location.end > b[0]:
                    seq.features[i].qualifiers['is_phage'] = 1
        except KeyError:
            pass
    sys.stdout.write(seq.format('genbank'))
