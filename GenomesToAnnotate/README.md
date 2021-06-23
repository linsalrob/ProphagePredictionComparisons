# Genomes To Annotate

This directory contains a selection of genomes that we are interested in manual curation. 

We have two sets of data for you to analyze:

- [less_hard](less_hard/) are genomes from the `Actinobacteria`, `Bacteroidota`, `Campylobacterota`, `Cyanobacteria`, `Desulfobacterota`, `Firmicutes`, `Proteobacteria`, or `Spirochaetota`. We have a lot of prophages already annotated in those genomes, and so we have a _better_ (but not complete) understanding of what they look like.
- [harder](harder/) is our diverse selection of bacterial genomes (i.e not those listed above). We are trying to diversify our phage selection, and so these are explicitly _not_ Proteobacteria, Fimicutes, or Bacteroidetes!

_Note:_ The genomes in the [less_hard](less_hard/) directory are a slightly different format: They are `gzip` compressed `gff3` format files. If you are using [Artemis](https://www.sanger.ac.uk/tool/artemis/) to annoate the genomes you do not need to uncompress them before opening them: Artemis will open them directly from the compressed format.

Please take a while and look thorugh these genomes. They have been automatically annotated with PhiSpy and we marked genes that have hits to pVOGs and other genes that might be phage genes. We have marked regions of the genome
that we think are phages.

## Getting started

If you are at [Rob's repository](https://github.com/linsalrob/ProphagePredictionComparisons/tree/ProphageAnnotations) click the ![fork me](../img/fork.png) button on the top right and copy this repository to your own GitHub. 

You don't need to clone the repository! You can download a file, manually annotate it, and upload it by dragging and dropping it on the web page. Add a little message, and you are good to go!


Once you have manually checked one of these genomes, please move it into the [AnnotatedGenomes](AnnotatedGenomes) folder, and then commit that to your repository. Whenever you are ready, please make a pull request to Rob's main repository, and 
we will grab the new genomes.

The main thing we would like is a `misc_feature` that delineates the end of the prophage genomes. Once we have that, we can munge the file into a suitable format.

Notes:
1. We have provided these as zip files for ease of use on Windows. Please commit the file in what ever format you like, we can handle it!
2. We recommend (still) using [Artemis](https://www.sanger.ac.uk/tool/artemis/) to annoate the genomes.
3. Pleae upload it to your repo and make a pull request:
   i. It allows us to credit you for your work
   ii. It allows us to keep track of new changes and integrate them
   iii. It is good practice!


