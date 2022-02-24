import sys
import gzip
from Bio import SeqIO

# use this script to get a table of the prophage coordinates from the genbank files
# e.g.:
# $ for i in `ls | grep gb.gz`; do python ../scripts/genbank2table.py $i; done > prophageCoords.tsv

orgName = {}
with gzip.open(sys.argv[1], 'rt') as fh:
    for line in fh:
        l = line.strip().split()
        if l[0] == 'DEFINITION':
            orgName['D'] = line.replace('DEFINITION','').strip()
        elif l[0] == 'SOURCE':
            orgName['S'] = line.replace('SOURCE', '').strip()
        elif l[0] == 'ORGANISM':
            orgName['O'] = line.replace('ORGANISM', '').strip()
            break

if len(orgName['O']) > 16:
    oName = orgName['O']
elif len(orgName['S']) > 16:
    oName = orgName['S']
else:
    oName = orgName['D']

fh = gzip.open(sys.argv[1], 'rt')
genbankfilesarestupid = SeqIO.parse(fh, 'genbank')
pp = 0
pgnum = 0
ptog = False
pstart = 0
pstop = 0
psid = ''
for seq in genbankfilesarestupid:
    psid = seq.id
    for i in seq.features:
        if i.type == 'CDS':
            try:
                i.qualifiers['is_phage'] == '1'
                pstop = str(i.location.end)
                pgnum += 1
                if ptog == False:
                    pp += 1
                    pstart = str(i.location.start)
                    ptog = True
            except KeyError:
                if ptog == True:
                    sys.stdout.write('\t'.join([sys.argv[1],oName,psid,str(pp),pstart,pstop,str(pgnum)]) + '\n')
                    ptog = False
                    pgnum = 0
                    pstart = 0
                    pstop = 0
    if int(pstart) > 0:
        sys.stdout.write('\t'.join([sys.argv[1],oName, psid, str(pp), pstart, pstop, str(pgnum)]) + '\n')
        ptog = False
        pgnum = 0
        pstart = 0
        pstop = 0

