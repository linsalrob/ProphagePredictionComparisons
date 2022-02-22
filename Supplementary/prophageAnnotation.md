# Manual curation of prophage regions in bacterial genomes

_Przemysław Decewicz (1)_, _Michael Roach (2)_

<small>1. Department of Environmental Microbiology and Biotechnology, Institute of Microbiology, Faculty of Biology, University of Warsaw, Warsaw, 02-096, Poland<br></small>
<small>2. Flinders Accelerator for Microbiome Exploration, Flinders University, 5042, SA, Australia</small>

## Looking for candidate regions 

Start by identifying proteins with sequence homology to known viral (phage) proteins and incorporate these annotations into your favourite genome browser (for instance, modifying CDS colors in your genbank file to highlight phage genes or as BED or GFF annotations). Inspect the genome and look for regions with the following: 

- Annotations including typical phage annotation words: phage, viral, virus, integrase, recombinase, transposase, repressor, terminase, capsid, coat, head, coral, head, neck, tail, sheath, baseplate, fiber, portal, morphogenesis, protease, holin, lysin, endolysin, ner, gemA, Mu, etc. 
- Tracks of hypothetical proteins 
- Tracks of genes transcribed in the same direction, often short or shorten than the surrounding ones, expect for Caudovirales phages where we may expect to see a few long CDSs of tail proteins 
- Local changes in GC% content over the potential prophage region 
- Presence of direct repeats surrounding potential prophage region (indicating potential attachment sites), which often occur in tRNA genes (these happen to colocalize next to tyrosine recombinase genes) but also in intergenic regions or within genes, eg. dusA 

Once a candidate region is observed manual inspection includes searching for the presence of 'the most conserved proteins' a prophage should have or we expect to find (these will be different for various Caudovirales families and Inoviruses), which includes a protein responsible for integration, lysis-lysogeny switch control, DNA packaging, major capsid protein, other structural protein-encoding genes, and lysis. These are found by homology searches against various databases, including NR (also limited to Viruses on taxonomy level) and CDD to identify domains that could suggest the viral origin of the protein.

In case we fail to identify any of these using NR and CDD (in practice it is either quick-blast limited to 50 hits or blastp limited to 50 hits with 1e-10 threshold), an additional search with the HHpred is performed (this one is very sensitive).

## Determining the range of prophage genome

In case when the direct repeats are not found, the prophage region ends are not clearly identifiable by the large change of GC%, presence of obvious bacterial host's metabolic genes or the potential region is Mu-like phage it is also possible to make the following searches: 

- Try to find nearly identical prophage region in a different bacterial genome in NT database and search for the presence of potential of attachment sites there, as well as check other options 
- Try to find a genome that does not have a prophage integrated, for this you can either: 
   - Take a prophage region with additional 5-10k regions surrounding the prophage candidate and search NT database - if the genome you're searching is deposited the first hit covers the whole region but the next ones may have a gap in the part where a prophage occurs - inspect that genome 
   - Take just the parts surrounding the candidate prophage region, join them into a single continuous sequence and search NT database, you may either get a perfect hit or with a gap/overlap depending on whether you left a piece of a prophage region or you added some bacterial genome part to prophage region 
   - Depending on the number of available genomes you can skip the first search and move to the second 

## When the prophage is incomplete 

It is better to indicate even a part of a prophage we are certain about as it distinguishes bacterial and phage genomes and allows for better training of classifiers.
Phage remnants are indicators of past infections and may increase the diversity of known phage-encoded proteins. 

## Feedback

These are some guidelines that we've written to describe the process we use to manually curate prophage regions. 
If you have any suggestions to improve this document we would love to hear from you!
Simply open an issue with your suggestions; or fork, edit, and create a pull request.