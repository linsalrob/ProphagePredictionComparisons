# Installation

Install conda, and then install snakemake (via conda is easiest).
Installation of all other software is handled by the snakemake pipelines. 
When running the pipelines, make sure you are using the flag `--use-conda`.

# Manual installation steps (none of the below steps should be necessary)

These instructions should enable you to set up your system to run the different comparisons.

For all the installations, we use [conda](https://docs.conda.io/en/latest/) to create virtual
environments, each one unique to the application. This allows us to sandbox out the different
python versions (e.g PhiSpy uses Python 3.8 and Seeker uses 3.7).

Most of these installation instructions come from the individual websites, so be sure to
check if they have been updated.

For each application we make a unique environment. You can list your current environments
with:

```
conda info -e
```

### Setting PERL5LIB

*Note:* Both [VirSorter](#VirSorter) and [ProphET](#prophet) suffer from a [known issue](https://stackoverflow.com/questions/58290190/how-to-fix-perl-from-anaconda-not-installing-bioperl-bailing-out-the-installat) 
with `conda` and `perl` that the libraries are not stored in exactly the right place. I did not do a deep dive into this because the solution is relatively simple, you need to modify `PERL5LIB`. Using
the solution above, you need to find the location of your conda environments (`conda info | grep envs`) and then find the PERL libraries in the appropriate environment. Finally, you need to modify the
`PERL5LIB` environment variable. In my instance it is to add this:

```bash
export PERL5LIB=$PERL5LIB:$HOME/anaconda3/envs/prophet/lib/perl5/site_perl/5.22.0/:$HOME/anaconda3/envs/virsorter/lib/perl5/site_perl/5.22.0/
```

Below, we show you how to set that for each environment, so you don't need to set it globally!


# [PhiSpy](https://github.com/linsalrob/PhiSpy)

```bash
conda create -n phispy -c bioconda phispy
```

# [VirSorter](https://github.com/simroux/VirSorter)

Please check the website above in case these instructions have changed. Please also see the note above about setting your `PERL5LIB` environment variable.

Install the databases

```bash
INSTALLDIR=~/virsorter
mkdir $INSTALLDIR
cd $INSTALLDIR
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
echo "dd12af7d13da0a85df0a9106e9346b45 virsorter-data-v2.tar.gz" > virsorter-data-v2.md5
md5sum -c virsorter-data-v2.md5
```

This should return the message: `virsorter-data-v2.tar.gz: OK`

```bash
tar -xvzf virsorter-data-v2.tar.gz
```

```bash
conda create --name virsorter -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond=0.9.14
conda install --name virsorter -c bioconda metagene_annotator
git clone https://github.com/simroux/VirSorter.git
cd VirSorter/Scripts
make clean
make
```


Use `conda info | grep envs` to get your conda environments directory, and then change to `virsorter/bin` inside that environment directory and set these two symbolic links.

```bash
ln -s $INSTALLDIR/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl
ln -s $INSTALLDIR/VirSorter/Scripts
```

Set the `PERL5LIB` environment variable. At this point **make sure you change the path to suit your environment**!

```
conda activate virsorter
conda env config vars set PERL5LIB=$PERL5LIB:$HOME/anaconda3/envs/virsorter/lib/perl5/site_perl/5.22.0/
conda deactivate
```

# [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/)

```bash
conda create -n checkv  -c conda-forge -c bioconda checkv
```

*Note*: According to [the README](https://bitbucket.org/berkeleylab/checkv/src/master/) CheckV has not been tested to predict prophages.

# [Seeker](https://github.com/gussow/seeker)

```bash
conda create --name seeker python=3.7 pip
conda activate seeker
pip install seeker
```

*Note*: According to [the README](https://github.com/gussow/seeker#python-library) Seeker is not suitable for predicting prophages.


# [ProphET](https://github.com/jaumlrc/ProphET/)

Please see the note above about setting your `PERL5LIB` environment variable.

```bash
conda create -n prophet -c bioconda blast-legacy emboss bedtools perl-bioperl perl-lwp-simple perl-gd
conda activate prophet
```

Set the `PERL5LIB` environment variable. At this point **make sure you change the path to suit your environment**!
We also deactivate/reactivate the environment to load it.

```
conda env config vars set PERL5LIB=$PERL5LIB:$HOME/anaconda3/envs/prophet/lib/perl5/site_perl/5.22.0/
conda deactivate
conda activate prophet
```

Continue with the install:
```
INSTALLDIR=~/prophet
mkdir $INSTALLDIR && cd $INSTALLDIR
git clone https://github.com/jaumlrc/ProphET.git
cd ProphET
./INSTALL.pl
```

*Note*: At this point, the installation currently (6/20/20) crashes. There are some [issues on GitHub](https://github.com/jaumlrc/ProphET/issues/30#issuecomment-647040938) that suggest it is because ProphET uses `http` and not `https` to download data from GenBank,
but it is not clear if that is the cause.


# [Phigaro](https://github.com/bobeobibo/phigaro)

```bash
INSTALLDIR=~/phigaro
mkdir $INSTALLDIR && cd $INSTALLDIR
conda create --name phigaro -c bioconda python=3.7 pip prodigal hmmer
conda activate phigaro
pip install phigaro
phigaro-setup --no-updatedb
```

*Notes*: 
- The download requests root access unless you provide the --no-updatedb option
- The download gets an undocumented file (which *appears* to contain HMMs) from a server in Russia
- The download uses uncompressed files of an undetermined size, and is very, very slow (about 1 kb every 2 minutes)



* [Phage_Finder](http://phage-finder.sourceforge.net/)

This is an old piece of software that is still available, so lets include it too!

```bash
INSTALLDIR=~/phage_finder 
mkdir $INSTALLDIR && cd $INSTALLDIR
conda create -n phage_finder -c bioconda blast-legacy hmmer mummer trnascan-se aragorn perl-math-round
wget https://downloads.sourceforge.net/project/phage-finder/phage_finder_v2.1/phage_finder_v2.1.tar.gz
echo "5f35122dced9438f5bf2798bddce156a  phage_finder_v2.1.tar.gz" > phage_finder_v2.1.md5
md5sum -c phage_finder_v2.1.md5
```

The md5sum command should say: `phage_finder_v2.1.tar.gz: OK`

```bash
tar xf phage_finder_v2.1.tar.gz
```

If you do not install in your `$HOME` directory, you will need to edit three files, `bin/HMM3_searches.sh`, `bin/phage_finder_v2.1.sh`, and `bin/Phage_Finder_v2.1.pl` to change
the value of `phome` to point to the appropriate installation location.

Activate the environment and set the PERL5LIB variable as above:

```
conda activate phage_finder
conda env config vars set PERL5LIB=$PERL5LIB:/home3/redwards/opt/phage_finder/phage_finder_v2.1/lib/
conda deactivate
```


