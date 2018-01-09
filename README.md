# Lahontan

Lahontan is a set of tools built around an object-oriented pipeline used in
qualifying and processing sets of paired-end next-gen RNA sequencing data. It
takes fastq files from a high-throughput sequencer and performs quality control
steps before performing a preliminary differential gene expression analysis.

Lahontan uses trimmomatic and fastqc to perform quality control; BLAST to check
the strandedness; and a choice of hisat2/featureCounts, hisat2/stringtie,
bowtie2, or kallisto to analyze and quantify the experimental data. Lahontan
also automatically produces some R scripts to begin a gene expression analysis
for certain modes. Lahontan is compatible with SLURM scheduled clusters for
parallel execution as well as standard desktop computers.

Table of Contents
=================

   * [Table of Contents](#table-of-contents)
   * [Installation](#installation)
   * [Built With](#built-with)
   * [Creating a Project and Running the Pipeline](#creating-a-project-and-running-the-pipeline)
      * [Requirements](#requirements)
      * [Setup](#setup)
      * [Quick Start](#quick-start)
         * [Example Execution](#example-execution)
   * [Pipeline Specifications](#pipeline-specifications)
      * [Project Stages](#project-stages)
      * [Project Structure](#project-structure)
      * [Usage: lahontan [options] [command] [args...]](#usage-runpipe-options-command-args)
         * [Available Commands](#available-commands)
         * [Command Options](#command-options)
            * [fcounts](#fcounts)
            * [string](#string)
            * [kall](#kall)
            * [bowtie2](#bowtie2)
            * [mj](#mj)
            * [mb](#mb)
            * [clean](#clean)
            * [fo](#fo)
            * [prepref](#prepref)
            * [gendef](#gendef)
            * [help](#help)
   * [Authors](#authors)

# Installation

To install, please see the [Wiki](https://github.com/UNR-CFB/rna-seq/wiki) for
detailed build instructions.

To summarize the [Wiki](https://github.com/UNR-CFB/rna-seq/wiki):

```
git clone --recursive https://github.com/UNR-CFB/lahontan.git
cd lahontan
./lib/autoSetup.sh
source ~/.bashrc
```
OR
```
singularity pull --name lahontan shub://UNR-CFB/lahontan
```

# Built With

* [Python](https://www.python.org/) (>=3.5)
* [Bash](https://www.gnu.org/software/bash/)
* [R](https://www.r-project.org/)
* [Kallisto](https://pachterlab.github.io/kallisto/)
* [Stringtie](https://github.com/gpertea/stringtie)
* [featureCounts](http://subread.sourceforge.net/)
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [hisat2](https://github.com/infphilo/hisat2)
* [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
* [samtools](https://github.com/samtools/samtools)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Seqtk](https://github.com/lh3/seqtk)
* [gffcompare](https://github.com/gpertea/gffcompare)
* [gffread](https://github.com/gpertea/gffread)

# Creating a Project and Running the Pipeline

## Requirements

* Paired-end Illumina fastq files. This means 2 files per sample. One for
  Read1, the other for Read2
* A reference genome for the organism in study. All chromosomes in a single
  FASTA format file
* A GTF-format file for containing the gene and exon coordinates for your
  reference genome
* A reference transcriptome (or some set of known transcripts) for your
  organism
* The Pipeline directory should be in your `$PATH` environmental variable. This
  can be accomplished by being in the main rna-seq directory and running:
```
$ cd scripts/Pipeline
$ echo 'export PATH=$PATH:'`pwd` >> $HOME/.bashrc
$ source $HOME/.bashrc
```

## Setup

It is important that the Input file that you feed into the pipeline follow the
format specified. You can use the `gendef` command to create a default input
file for you, which at that point only requires you to fill in the Project,
Reference, and Original variable directories. See the [`gendef` documentation](#prepref)
for more help. An example Input file can be seen below:

```
ExampleInputFile
-------------------------------------------------------------------------------
[locations]
Project = /home/user/Projects/ProjectName
Reference = /home/user/References/ProjectName-Reference
Original = /home/user/Datas/ProjectName-Data
```

Often times your data will be organized in a way that isn't accepted by
lahontan.  In order to continue, your `Original` directory must be arranged such
that each sample read is alphanumerically sorted in a single directory (Note:
sorted by sample and by read)

Example Original Directory:
```
/home/user/Data/ProjectName-Data/
    SpeciesName.treament1_group1.read1.fq.gz
    SpeciesName.treament1_group1.read2.fq.gz
    SpeciesName.treament1_group2.read1.fq.gz
    SpeciesName.treament1_group2.read2.fq.gz
    SpeciesName.treament2_group1.read1.fq.gz
    SpeciesName.treament2_group1.read2.fq.gz
    SpeciesName.treament2_group2.read1.fq.gz
    SpeciesName.treament2_group2.read2.fq.gz
```

Also, in order to limit errors, your filenames must only contain: letters,
numbers, periods, hyphens, and underscores.

Your reference filenames must also follow a few specific rules. Along with only
using letters, numbers, periods, hyphens, and underscores:
* Your GTF-format file must end with ".gtf" and it must be the only file within
  your reference directory to do so
* Your reference transcriptome must contain ".cdna." anywhere within the name
  except for the end. It must also be the only file within your reference
  directory to do so
* Your reference genome must contain ".dna." anywhere within the name except
  for the extension.  It must also be the only file within your reference
  directory to do so

Example Reference Directory:
```
/home/user/Reference/ProjectName-Reference/
    SpeciesName.cdna.all.fa.gz
    SpeciesName.dna.toplevel.fa.gz
    SpeciesName.gtf
```

## Quick Start

1. Create an Input file with `lahontan gendef` for your desired pipeline using
   the `-m` option. Then fill in the "[locations]" section.
2. Make sure data and references are following format specified in Setup
3. Prepare reference data with `lahontan prepref`. If you are using the `kall`
   or `bowtie2` pipeline, be sure to add the `--kallisto` or `--bowtie2` flag
   to make sure the appropriate indices get created.
4. Run `lahontan mj` to make a metadata file
5. Run `lahontan mb` to make a blacklist file
6. Decide if you will be running the pipeline with SLURM or not
    * If yes, then run `lahontan fo` with proper arguments to make optimal path
      file
    * If no, then you can skip this step
    * The following steps will be denoted with a (s) if the step is for SLURM or
      a (w) if the step is without SLURM
7. (s) Make a batch file by adding the `--makebatch` option to any of the three
   pipeline commands. You should also add the `--batchjson` option with the
file you made with `lahontan fo`; the `--jsonfile` option with the file you made
with  `lahontan mj`; the `--use-blacklist` option with the file you made with
`lahontan mb`; and the `--use-reference` because the reference data has been
prepared. So an example command to make a featureCounts batch file may look like:
```
$ lahontan fcounts --makebatch 48,48,32 --batchjson /path/to/OptimalPath.dat \
--jsonfile /path/to/Metadata.json --use-blacklist /path/to/CustomBlacklist.fa \
--use-reference /path/to/INPUT
```
7. (w) You can execute the entire pipeline with one command, but it's important
   that the `--jsonfile`, `--use-blacklist`, and `--use-reference` options are
added to the command. An example command may look like:
```
$ lahontan fcounts --jsonfile /path/to/Metadata.json \
--use-blacklist /path/to/CustomBlacklist.fa --use-reference --maxcpu 24 \
/path/to/INPUT
```
8. (s) Execute the batch file with `sbatch pipeBatch`

Note: It may be a good idea to hold your input file, metadata file, blacklist
file, batch file, and optimal path file in the same directory for easy tracking

### Example Execution

Assuming that the data and reference files have been appropriately organized,
and are located at `/home/user/Data/MyDataFolder` and
`/home/user/References/MyReferenceFolder`, run the following to obtain a
default input file:
```
$ lahontan gendef -m fcounts -f MyFcountsProject.ini
```
At this point complete the "[locations]" section in `MyFcountsProject.ini`. So
`MyFcountsProject.ini` at minimum contains:
```
MyFcountsProject.ini
-------------------------------------------------------------------------------
[locations]
Project = /home/user/Projects/MyProjectFolder
Reference = /home/user/References/MyReferenceFolder
Original = /home/user/Datas/MyDataFolder
```
Then, execute the following:
```
$ lahontan prepref -r /home/user/References/MyReferenceFolder
$ lahontan mj -j MyMetadata.json
$ lahontan mb -t MyBlacklist.fa
$ lahontan fcounts --use-reference --use-blacklist MyBlacklist.fa \
--jsonfile MyMetadata.json MyFcountsProject.ini
```

# Pipeline Specifications

## Project Stages

**Stage 1**

Will create the [project structure](#project-structure) and will link the
original data files, the reference files, and a Metadata file if it has been
prepared. This stage is identical for all pipelines.

**Stage 2**

Will run a quality check on the reference data and preprocess the reference
data. This entails building a BLAST database, finding the splice sites, finding
the exons, building a hisat2 database, and building a genome index with
samtools. This stage can be identically run with `lahontan prepref` outside of a
Project structure and before running a pipeline. Stringtie and featureCounts
will have an identical Stage 2, but Kallisto requires the construction of its
own index. If preparing the references beforehand for Kallisto, make sure to
add the `--kallisto` flag to the `lahontan prepref` command.

**Stage 3**

The first step for all pipelines in Stage 3 is to run a quality check on the
data. This includes using fastQC to create statistics on the quality of the
data as well as to search for overrepresented sequences. Trimming of the reads
is done by trimmomatic which trims poor reads as well as any reads in a
blacklist file. A custom blacklist file can be created by `lahontan mb`, which
can be fed to a pipeline with the `--use-blacklist` option. After trimming, an
additional round of fastQC will be performed to check the results of the
trimming. These quality control steps are identical for all pipelines. Once the
quality control is finished, the strandedness of the reads is checked with
BLAST.

After quality control is done, depending on which pipeline is being ran, the
next steps will vary. For featureCounts, hisat2 will be used to align the reads
to the references, then featureCounts will quantify the gene expression. The
pipeline will then gather the data to create tables which can be used in
further analysis.

If Stringtie is being used, then hisat2 will be used to align the reads to the
references. After the alignment, Stringtie will build a new reference transcripts
file which will then be used in estimating the transcript abundances.

If Kallisto is being used, then Kallisto will perform its probabilistic
quantification on the raw data using the index it built.

**Stage 4**

In Stage 4, the resulting data from whichever pipeline was ran is tabulated and
prepared for statistical analysis with some automatically generated R scripts.
The automatically generated R scripts are built with deseq2 and edgeR for
featureCounts, ballgown for Stringtie, and sleuth for Kallisto.

**Stage 5**

In Stage 5, the automatically generated R scripts from Stage 4 are executed.

## Project Structure

Stage 1 of the Pipeline will create a consistent directory structure between
Projects that includes a *Data*, *Original*, *Postprocessing*, and *Reference*
directory which looks like:
```
ProjectName/
    Data/
        sample_01/
            sample_01-read_1
            sample_01-read_2
        sample_02/
            sample_02-read_1
            sample_02-read_2
        sample_03/
            sample_03-read_1
            sample_03-read_2
        .
        .
        .
        (Other sample directories)
    Original/
        sample_01-read_1
        sample_01-read_2
        sample_02-read_1
        sample_02-read_2
        sample_03-read_1
        sample_03-read_2
        .
        .
        .
        (Other sample files)
    Postprocessing/
        (Metadata.json if already created, which I recommend you do)
    Reference/
        Genome.fa
        cDNA.fa
        .gtf
        (Other files if reference files have already been processed and prepared)
```
* *Data* = directory where Stage 3 will propogate. Contains a directory for each
respective sample which contains its respective raw data files
* *Original* = directory containing symbolic links pointing to original raw data files
* *Postprocessing* = directory where differential expression analysis and other analyses will be saved to
* *Reference* = directory containing symbolic links pointing to original reference files

## Usage: lahontan [options] [command] [args...]

The pipeline is designed as a series of scripts built around a central
function, `lahontan`. `lahontan` functions when it is fed a command or argument.
If it is fed a command, that activates the options available for that command,
and can subsequently be added after the command. The available commands which
can be given to `lahontan` are shown below.

### Available Commands

```
fcounts     For running featureCounts pipeline
string      For running Stringtie pipeline
kall        For running kallisto pipeline
bowtie2     For running Bowtie2 pipeline
mj          Provoke questionnaire to make a Metadata file
mb          Create trimmomatic blacklist
clean       Clean any project directories
fo          Finds optimal execution path for batch execution
prepref     Pre-processes reference data
gendef      Generate a default input file
```

For more information about any of the commands, you can add the `-h` option
to see the available options and any miscellaneous information. You can also give
the command as an argument to the `help` command.

For example:
```
$ lahontan fcounts -h
$ lahontan help fcounts
```
will show information available about the `fcounts` command. `lahontan` also has
various arguments which can be used to perform an action or change the behavior
of the pipeline. Shown below are the available options for `lahontan`.

**Options:**

```
-h, --help
    Show this screen and exit
--version
    Show version and exit
--noconfirm
    Ignore all user prompts except JSON file creation
```

The syntax for the options displayed is based on the conventional CLI
interface.  The list of options describe whether an option has short/long forms
(`-h, --help`), whether an option has an argument, and whether the argument has
a default value. If an option has both a short and long form, then either form
may be used for the same functionality. The default values will be denoted in
the description of the option by, `[default: foo]`.

**More information about options:**

--version
> Will show the version of `lahontan` currently installed

--noconfirm
> Used to avoid any interactive prompts that may be shown which would require
user intervention

### Command Options

The commands available to `lahontan` will be described in more detail below.

#### fcounts
**Usage: lahontan fcounts [options] INPUTFILE**

When `lahontan` is given the `fcounts` command, the pipeline is given the
information to use the featureCounts tools and pathway for analysis. If given
the fcounts command, then a mandatory argument `INPUTFILE` is required. It also
has various optional arguments which can be given to change the behavior of the
pipeline. The `INPUTFILE` argument expects the path to a file that defines
`Project`, `Reference`, and `Original` variables. The `Project` variable should
define the path to a directory where the pipeline will save and execute its
analyses. At first execution, the `Project` directory should be empty, as the
pipeline will populate the directory into a consistent project structure.  The
`Reference` variable should define the path to a directory where the reference
materials are stored; at minimum it should contain a gene transcripts file, a
cDNA file, and a reference genome. The `Original` variable should define the
path to a directory that contains the paired-end sequencing data.

An example of an `INPUTFILE`:
```
ExampleInputFile
-------------------------------------------------------------------------------
[locations]
Project = /home/user/ProjectName
Reference = /home/user/Reference/ProjectName-Reference
Original = /home/user/Data/ProjectName-Data
```

**Options**

```
-h, --help
    Show this screen and exit
-e <stage>, --execute <stage>
    Comma-separated list of stages to be executed.
    Possible stages include:
        1: Creating Project Structure
        2: Preparing Reference Data
        3: Running actual Pipeline
        4: Preparing for R analysis
        5: Running R analysis
        A: (1,2,3,4,5); A=all i.e. runs entire pipeline
    [default: A]
-r <integer>, --runsample <integer>
    Runs Stage 3 of the pipeline on the sample specified
    by the integer
-j <jsonFile>, --jsonfile <jsonFile>
    Ignores JSON Metadata file creation and uses specified
    path to JSON Metadata
--maxcpu <CPUs>
    Limits number of CPUs used by Pipeline. Default is to
    use all available CPUs
--makebatch <cluster>
    Makes batch file to be used with slurm. The argument
    it takes is a comma-separated list of CPUs on each
    node in your cluster
--batchjson <pathtoJSON>
    Uses json file already created to make batch file
--edger
    Runs edgeR analysis only. Default is to run both
--deseq
    Runs DESeq2 analysis only. Default is to run both
--noconfirm
    Ignore all user prompts except JSON file creation
--use-blacklist <blacklist>
    Trimmomatic blacklist used for quality control
    [default: $RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa]
--use-reference
    Use Reference data that has already been prepared.
    Put path to already prepared reference data in
    INPUT file
    Note: Need to run Stage 1 with this argument or add to
    "--makebatch" argument
```

**More information about options:**

-e *\<stage\>*, --execute *\<stage\>*
> Used to specify stage of Pipeline to be ran
> *\<stage\>* is a comma-separated list of stages to be executed.
> Possible stages include:
>> 1 : Creating Project Structure; this includes creating Project directories
>> and making symbolic links to reference and raw data
>> 2 : Preparing Reference Data; this includes running Quality Control on the
>> reference data then setting up hisat2 databases
>> 3 : Running actual Pipeline; this includes quality control, alignment,
>> compression, and feature counts
>> 4 : Preparing for R analysis; this includes creating a 'nice' count file and
>> developing automated R scripts from `Metadata.json`
>> 5 : Running R analysis; this involves running DESeq2 and edgeR differential
>> expression analyses
>> A : (1,2,3,4,5) i.e. runs entire pipeline

> Note: default is to just run `A`

-r *\<sampleNumber\>*, --runsample *\<sampleNumber\>*
> Used to execute through Stage 3 on the sample specified by *\<sampleNumber\>*

-j *\<jsonFile\>*, --jsonfile *\<jsonFile\>*
> Ignores JSON Metadata file creation and uses specified path to already
> created JSON Metadata

--maxcpu *\<CPUs\>*
> Used to limit number of CPUs used by Pipeline
> Note: default is to use all available CPUs

--makebatch *\<cluster\>*
> Used to make batch file for use with slurm
> *\<cluster\>* is a comma-separated list of CPUs on each node in your cluster

--batchjson *\<pathtoJSON\>*
> If an optimal path JSON file has been created using `lahontan fo` or another
> method, use option to skip the calculation of the optimal CPU usage path
> Note: Option is depended on `--makebatch` i.e. must be added with
> `--makebatch` or not at all

--edger
> Runs edgeR analysis only
> Note: default is to run both edgeR and deseq2 analyses

--deseq
> Runs deseq2 analysis only
> Note: default is to run both edgeR and deseq2 analyses

--noconfirm
> Used for non-interactivity
> Ignores all user prompts except JSON file creation unless given `--jsonfile` option

--use-blacklist *\<blacklist\>*
> Trimmomatic blacklist used for quality control. A new blacklist can be made with
> `lahontan mb`
> Note: default is to use `$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa` where
> `$RNASEQDIR` is the build directory specified in the installation instructions

--use-reference
> If reference data has already been prepared, use option to skip Stage 2 and
> use prepared reference data. It is important that this option be used for
> every stage if the reference data is already prepared.
> Note: If option is to be used, must be ran from Stage 1

#### string
**Usage: lahontan string [options] INPUTFILE**

Similar to the `fcounts` command, `string` gives the pipeline the information
it needs so that it uses the Stringtie tools and pathway. The options to the
`string` command are similar to those of `fcounts`. The only difference in
execution are Stages 3, 4, and 5. In stage 3, instead of running the analyses
in one step, Stringtie takes three phases to do so. The first phase, `a`, is
where transcripts are assembled for each sample in the experiment. The second
phase, `b`, is where a new gene transcript file is created from each sample's
transcripts. The third step, `c`, is where Stringtie estimates the transcript
abundances for each sample using the newly assembled gene transcript file.
Stage 3 of the pipeline can be ran altogether by using the `--phase abc`
option. In Stage 4, the pipeline prepares the results to be analyzed by
ballgown in Stage 5.

**Options**

```
-h, --help
    Show this screen and exit
-e <stage>, --execute <stage>
    Comma-separated list of stages to be executed.
    Possible stages include:
        1: Creating Project Structure
        2: Preparing Reference Data
        3: Running actual Pipeline
        4: Preparing for R analysis
        5: Running R analysis
        A: (1,2,3,4,5); A=all i.e. runs entire pipeline
    [default: A]
-r <integer>, --runsample <integer>
    Runs Stage 3 of the pipeline on the sample specified
    by the integer
-p <phase>, --phase <phase>
    Use stringtie tools to replace featureCounts. phase can
    be any of: "a","b","c","ab","bc","abc"
    Option also used to specify stringtie postprocessing
    options in which case use: "--execute 4 --stringtie abc"
    Note: If you will be running phase b, you cannot specify
    a sample to run with --runsample
    [default: abc]
-j <jsonFile>, --jsonfile <jsonFile>
    Ignores JSON Metadata file creation and uses specified
    path to JSON Metadata
--maxcpu <CPUs>
    Limits number of CPUs used by Pipeline. Default is to
    use all available CPUs
--noconfirm
    Ignore all user prompts except JSON file creation
--use-blacklist <blacklist>
    Trimmomatic blacklist used for quality control
    [default: $RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa]
--use-reference
    Use Reference data that has already been prepared.
    Put path to already prepared reference data in
    INPUT file
    Note: Need to run Stage 1 with this argument or add to
    "--makebatch" argument
--makebatch <cluster>
    Makes batch file to be used with slurm. The argument
    it takes is a comma-separated list of CPUs on each
    node in your cluster
--batchjson <pathtoJSON>
    Uses json file already created to make batch file
```

**More information about options:**

The `string` command only adds the `--phase` option from the `fcounts` command.
If you wish to see more information, see the [fcounts](#fcounts) documentation.

-p *\<phase>\*, --phase *\<phase>\*
> Modifies pipeline to incorporate stringtie and ballgown. Stringtie replaces
> featureCounts.
> Possible phases can be any of: "a","b","c","ab","bc","abc"
>> a : Phase a refers to mapping sample reads to reference genome, converting
>> mapped reads to .bam format, and then assembling transcripts for each sample
>> b : Phase b refers to merging the transcripts from each sample into one
>> merged .gtf file
>> c : Phase c refers to estimating the transcript abundances for each sample
> Note: If you will be running phase b, you cannot specify a sample to run with --runsample

#### kall
**Usage: lahontan kall [options] INPUTFILE**

Similar to the `fcounts` and `string` commands, `kall` gives the pipeline the
information it needs so that it uses the Kallisto tools and pathway. The
options to the `kall` command are similar to those of `fcounts`. The only
difference is that Kallisto is used instead of hisat2/featureCounts/Stringtie
for quantifying gene expression.

**Options**

```
-h, --help
    Show this screen and exit
-e <stage>, --execute <stage>
    Comma-separated list of stages to be executed.
    Possible stages include:
        1: Creating Project Structure
        2: Preparing Reference Data
        3: Running actual Pipeline
        4: Preparing for R analysis
        5: Running R analysis
        A: (1,2,3,4,5); A=all i.e. runs entire pipeline
    [default: A]
-r <integer>, --runsample <integer>
    Runs Stage 3 of the pipeline on the sample specified
    by the integer
-j <jsonFile>, --jsonfile <jsonFile>
    Ignores JSON Metadata file creation and uses specified
    path to JSON Metadata
--maxcpu <CPUs>
    Limits number of CPUs used by Pipeline. Default is to
    use all available CPUs
--noconfirm
    Ignore all user prompts except JSON file creation
--use-blacklist <blacklist>
    Trimmomatic blacklist used for quality control
    [default: $RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa]
--use-reference
    Use Reference data that has already been prepared.
    Put path to already prepared reference data in
    INPUT file
    Note: Need to run Stage 1 with this argument or add to
    "--makebatch" argument
--makebatch <cluster>
    Makes batch file to be used with slurm. The argument
    it takes is a comma-separated list of CPUs on each
    node in your cluster
--batchjson <pathtoJSON>
    Uses json file already created to make batch file
```

**More information about options:**

The `kall` command does not add any options compared to the `fcounts` command.
If you wish to see more information, see the [fcounts](#fcounts) documentation.

#### bowtie2
**Usage: lahontan bowtie2 [options] INPUTFILE**

Similar to the `fcounts`,`string`, and `kall` commands, `bowtie2` gives the
pipeline the information it needs so that it uses the bowtie2 aligner. The
options to the `bowtie2` command are similar to those of `kall`.

**Options**

```
-h, --help
    Show this screen and exit
-e <stage>, --execute <stage>
    Comma-separated list of stages to be executed.
    Possible stages include:
        1: Creating Project Structure
        2: Preparing Reference Data
        3: Running actual Pipeline
        4: Preparing for R analysis
        5: Running R analysis
        A: (1,2,3,4,5); A=all i.e. runs entire pipeline
    [default: A]
-r <integer>, --runsample <integer>
    Runs Stage 3 of the pipeline on the sample specified
    by the integer
-j <jsonFile>, --jsonfile <jsonFile>
    Ignores JSON Metadata file creation and uses specified
    path to JSON Metadata
--maxcpu <CPUs>
    Limits number of CPUs used by Pipeline. Default is to
    use all available CPUs
--noconfirm
    Ignore all user prompts except JSON file creation
--use-blacklist <blacklist>
    Trimmomatic blacklist used for quality control
    [default: $RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa]
--use-reference
    Use Reference data that has already been prepared.
    Put path to already prepared reference data in
    INPUT file
    Note: Need to run Stage 1 with this argument or add to
    "--makebatch" argument
--makebatch <cluster>
    Makes batch file to be used with slurm. The argument
    it takes is a comma-separated list of CPUs on each
    node in your cluster
--batchjson <pathtoJSON>
    Uses json file already created to make batch file
```

**More information about options:**

In the bowtie2 pipeline, the metadata file for into the `--jsonfile` option
serves only as a reference for the user and offers no use in any of the
`bowtie2` tools.

#### mj
**Usage: lahontan mj [options]**

Due to the many different types of experiments which can be run, it is is
difficult to automate a method for gathering information about the treatments
and groups applied in an experiment. The treatments and groups of an experiment
are important in being able to develop any conclusions with the data. The
`lahontan mj` command provides a questionnaire that requires user feedback about
the experiment. The important point is to provide enough information about the
features of the experiment so a sample can be identified by those features.

Caution: It is important to check that the order with which the samples are
linked in their respective directory is the same as the features you describe
it with.

For example, a typical experiment might look like:
```
$ lahontan mj
What is the name of the project? ProjectName
How many samples are there? 4
How many features does each sample have? 2
What is the name of feature #1? treatment
What is the name of feature #2? group
What is the name of the main feature? treatment
What is the treatment feature for sample #1? X
What is the group feature for sample #1? GroupA
What is the treatment feature for sample #2? Y
What is the group feature for sample #2? GroupA
What is the treatment feature for sample #3? X
What is the group feature for sample #3? GroupB
What is the treatment feature for sample #4? Y
What is the group feature for sample #4? GroupB
Done making Metadata.json
```
In this way, you can distinguish any sample by specifying features #1 and #2.
If there were multiple runs with the same treatment and group, you may choose
to add a "run" feature. This questionnaire creates a JSON file by default
titled `Metadata.json`:
```
Metadata.json
--------------------------------------------------------
{
    "FeatureNames": [
        "group",
        "treatment"
    ],
    "MainFeature": "treatment",
    "NumberofFeatures": 2,
    "NumberofSamples": 4,
    "ProjectName": "ProjectName",
    "Samples": {
        "sample_01": {
            "Features": {
                "group": "GroupA",
                "treatment": "X"
            }
        },
        "sample_02": {
            "Features": {
                "group": "GroupA",
                "treatment": "Y"
            }
        },
        "sample_03": {
            "Features": {
                "group": "GroupB",
                "treatment": "X"
            }
        },
        "sample_04": {
            "Features": {
                "group": "GroupB",
                "treatment": "Y"
            }
        }
    }
}
```

**Options**

```
-h, --help
    Show this screen
-j <jsonfile>, --jsonfile <jsonfile>
    Optional name of JSON file to be saved to [default: Metadata.json]
```

**More information about options:**

-j *\<jsonfile\>*, --jsonfile *\<jsonfile\>*
> Optional name of JSON file to be saved to
> Note: default is to save the file to `Metadata.json` in the current directory

#### mb
**Usage: lahontan mb [options]**

This command is used to create a blacklist of reads that qualify as noise and
which should be trimmed from the data by trimmomatic.

As of trimmmomatic version 0.36, the default blacklist file is:
```
$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE.fa
-------------------------------------------------------------------------------
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

The blacklist file created by `lahontan mb`, as of lahontan version 1.1 is:
```
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>ABISolid3AdapterB
CCTATCCCCTGTGTGCCTTGGCAGTCTCAGCCTCTCTATGGGCAGTCGGT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA
>ConsecutiveA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>ConsecutiveT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ConsecutiveC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>ConsecutiveG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
```

**Options**

```
-h, --help
    Show this screen and exit
-t <newadapterfile>, --tofile <newadapterfile>
    Where to save new blacklist [default: ./CustomBlacklist.fa]
```

**More information about options:**

-t *\<newadapterfile\>*, --tofile *\<newadapterfile\>*
> Optional name of file to be saved to
> Note: default is to save the file to `CustomBlacklist.fa` in the current directory

#### clean
**Usage: lahontan clean [options] INPUTFILE**

This command is used to clean/reset the populated project directories. The
default behavior is to clean out the entire Project structure to the point
right after Stage 1.

**Options**

```
-h, --help
    Show this screen and exit
-c <placeToClean>, --clean <placeToClean>
    Cleans <placeToClean>; Possible places include:
        Reference
        Data
        Postprocessing
        All
--sampleclean <sampleName>
    Similar to --clean; but instead just cleans a
    single sample directory, <sampleName>
```

**Mandatory Arguments**

INPUTFILE
> The path to a file that contains the variables `Project`, `Reference`, and
> `Original` defined by their respective paths

**More information about options:**

-c *\<placeToClean\>*, --clean *\<placeToClean\>*
> Used to clean/reset *\<placeToClean\>*
> Possible places include:
>> Reference = removes all files except Genome, cDNA, and gtf in Reference
>> directory
>> Data = removes all files except raw data file links in Data directory
>> Postprocessing = removes all files except `Metadata.json` if it exists in
>> Postprocessing directory
>> All = resets entire Project. This is equivalent to running stage 1 on fresh
>> data

--sampleclean *\<sampleName\>*
> Used to clean/reset *\<sampleName\>*
> Similar to `--clean Data` but instead of cleaning all samples just cleans a
> single sample directory
> Note: *\<sampleName\>* will be of the form `sample_#` where #
> refers to specific sample number


#### fo
**Usage: lahontan fo [options] NUMSAMPLES CLUSTER**

It is possible to run the pipeline in a SLURM cluster by creating a batch file
that contains various `srun` commands given to the SLURM scheduler. However,
the automated batch file maker needs information on how best to run the samples
in order to minimize the waste of CPU resources. The `lahontan fo` command
offers the ability to create a file that contains information on how best to
run the samples.


**Options**

```
-h, --help
    Show this screen and exit
-t <filename, --tofile <filename>
    Name of file to be saved to
    [default: ./OptimalPath.dat]
--maxcpu <CPUs>
    Limits number of CPUs used in calculating optimal path
    Default is to use all available CPUs
-c, --customize
    If the calculation is too slow, can use this option
    to specify the path yourself
```

**Mandatory Arguments**

NUMSAMPLES
> The number of samples in your experiment.

CLUSTER
> A comma separated list of the number of CPUs on each node in your cluster

**More information about options:**

-t <filename, --tofile *\<filename\>*
> Optional name of file to be saved to
> Note: default behavior is to save to OptimalPath.dat in the current directory

--maxcpu <CPUs>
> Limits number of CPUs used in calculating optimal path
> Note: Default is to use all available CPUs

-c, --customize
> If the calculation is too slow, can use this option to specify the path
> yourself

#### prepref
**Usage: lahontan prepref [options] (-r REFERENCEDIR | -i INPUT)**

This command is equivalent to Stage 2 of the pipeline, however the difference
is that this command is occurring outside of a structured Project environment.
It is useful to prepare the reference information before running the pipeline
as there may be unforeseen difficulties with the reference information which
can't be caught by the pipeline.

**Options**

```
-h, --help
    Show this screen and exit
-r <referenceDir>
    One of the mandatory arguments i.e. must specify either
    r OR i options.
    Option value is a path to directory that contains a GTF,
    cDNA, and reference genome
-i <input>
    One of the mandatory arguments i.e. must specify either
    r OR i options
    Option value is a path to input file that must contain
    at least a path to reference directory
    Note: a default input file can be generated with:
        lahontan gendef --mode prepref
-q, --qualitycheck
    Only run quality check on references. Default behavior
    is to run both quality control and preprocessing
-p, --preprocess
    Only run preprocessing on references. Default behavior
    is to run both quality control and preprocessing
-k, --kallisto
    Additionally build kallisto index
    Note: index required if kallisto pipeline to be used
--onlykallisto
    Build only the kallisto index
-b, --bowtie2
    Additionally build Bowtie2 index
    Note: index required if Bowtie2 pipeline to be used
--onlybowtie2
    Build only the Bowtie2 index
--maxcpu <CPUs>
    Limit the number of CPU that get used to preprocess
    reference data.
    Note: default is to use all available
```

**Mandatory Arguments**

(-r REFERENCEDIR | -i INPUT)
> For the `prepref` command, you must use either the `-r` option with the
> path to a directory that contains a GTF, cDNA, and reference genome OR
> you must use the `-i` option with a manifest file that includes at least
> the path to a reference directory that contains a GTF, cDNA, and reference
> genome

An example of a minimal input file for the `prepref -i` option:
```
ExampleInputFile
-------------------------------------------------------------------------------
[locations]
Reference = /home/user/Reference/ProjectName-Reference
```
Note the "[locations]" header and the "Reference" variable. These are required
in every `prepref` input file.

**More information about options:**

-q, --qualitycheck
> Check the quality of the references by comparing the chromosome names in the
> genome and gtf, comparing the gene names in the cDNA and gtf, and comparing
> the number of unique gene names in the cDNA and gtf.
> This option is given if you only want to check the quality without actually
> preprocessing the references

-p, --preprocess
> Preprocesses the references by building a BLAST database, finding the splice
> sites, finding the exons, building a hisat2 database, and building a genome
> index with samtools
> This option is given if you only want to preprocess the references without
> actually running a quality check

-k, --kallisto
> Build an index for kallisto analysis

-b, --bowtie2
> Build an index for bowtie2 analysis

--maxcpu <CPUs>
> Limit the number of CPU that get used to preprocess reference data.
> Note: default is to use all available

#### gendef
**Usage: lahontan gendef [options]**

**Options**

```
-h, --help
    Show this screen and exit
-f <filename>, --filename <filename>
    Name of input file to create
    [default: default_input.ini]
-m <mode>, --mode <mode>
    Type of execution mode, can be:
        fcounts
        string
        kall
        bowtie2
        prepref
    [default: fcounts]
```

The `gendef` command is used to generate an input file with all default
parameters filled in. An input file is mandatory for the `fcounts`, `string`,
`kall`, and `bowtie2` commands. It is optional for the `prepref` command. The
input files are organized in an INI format. For more information about file
format and syntax, see the [python
documentation](https://docs.python.org/3/library/configparser.html#supported-ini-file-structure)
. To summarize, files consist of sections, each with a section header in square
brackets, and key-value pairs defined underneath (separated by an "="). To
leave a comment in a file, you can precede a line with a "#" or ";".

For example, here is an input file that is generated with
```
$ lahontan gendef -m fcounts
```

```
default_input.ini
-------------------------------------------------------------------------------
[locations]
Project = /home/alberton/rna-seq/scripts/Pipeline/Test
Reference = /home/alberton/rna-seq/scripts/Pipeline/Test
Original = /home/alberton/rna-seq/scripts/Pipeline/Test

[mode]
fcounts = True
kall = False
string = False
bowtie2 = False

[fcounts]
--execute = A
--runsample = None
--jsonfile = None
--maxcpu = None
--noconfirm = False
--makebatch = None
--batchjson = None
--edger = False
--deseq = False
--use-blacklist = /home/alberton/rna-seq/build/Trimmomatic/adapters/TruSeq3-PE.fa
--use-reference = False

[EXECUTABLES]
Rscript = /usr/local/bin/Rscript
java = /usr/bin/java
rnaseqdir = /home/alberton/rna-seq/build
fastqc = ${rnaseqdir}/fastqc
extract_exons.py = ${rnaseqdir}/hisat2/extract_exons.py
extract_splice_sites.py = ${rnaseqdir}/hisat2/extract_splice_sites.py
hisat2 = ${rnaseqdir}/hisat2/hisat2
hisat2-build = ${rnaseqdir}/hisat2/hisat2-build
samtools = ${rnaseqdir}/samtools/bin/samtools
makeblastdb = ${rnaseqdir}/ncbi-blast/bin/makeblastdb
blastn = ${rnaseqdir}/ncbi-blast/bin/blastn
seqtk = ${rnaseqdir}/seqtk
gffcompare = ${rnaseqdir}/gffcompare
stranded_classifier.py = ${rnaseqdir}/stranded_classifier.py
featureCounts = ${rnaseqdir}/subread/bin/featureCounts
stringtie = ${rnaseqdir}/stringtie
kallisto = ${rnaseqdir}/kallisto
bowtie2 = ${rnaseqdir}/bowtie2
bowtie2-build = ${rnaseqdir}/bowtie2-build

[fcounts/main]
runQCheck = True
runTrimmomatic = True
runSeqtk = True
runBlastn = True
runHisat = True
runCompression = True
runFeatureCounts = True
getNiceColumns = True
getAlignedColumn = True

[fcounts/runHisat]
hisat2 = ${EXECUTABLES:hisat2}
-k = None
-p = None
--rna-strandedness = None
--dta = None
--phred = None
--known-splicesite-infile = None
-x = None
-1 = None
-2 = None
-S = None
other = None

[fcounts/runCompression]
samtools = ${EXECUTABLES:samtools}
-b = None
-T = None
-@ = None
-o = None
in = None
other = None

[fcounts/runFeatureCounts]
featureCounts = ${EXECUTABLES:featureCounts}
-T = None
-s = None
-p = None
-C = None
--primary = None
--ignoreDup = None
--largestOverlap = None
-t = None
-g = None
-a = None
-o = None
in = None
other = None

[runQCheck]
fastqc = ${EXECUTABLES:fastqc}
-t = None
-o = None
read1 = None
read2 = None
other = None

[runTrimmomatic]
java = ${EXECUTABLES:java}
-jar = None
-threads = None
-phred = None
read1 = None
read2 = None
read1pout = None
read2pout = None
read1uout = None
read2uout = None
blacklist = None
ILLUMINACLIP = None
LEADING = None
TRAILING = None
SLIDINGWINDOW = None
MINLEN = None
other = None

[runSeqtk]
seqtk = ${EXECUTABLES:seqtk}
-s = None
read1 = None
read2 = None
read1out = None
read2out = None
samples = None
-A = None
other = None

[runBlastn]
blastn = ${EXECUTABLES:blastn}
read1 = None
read2 = None
out1 = None
out2 = None
-db = None
-task = None
-outfmt = None
-max_target_seqs = None
-num_threads = None
other = None
```

This input file contains the default behavior of the `fcounts` pipeline.The
order of the sections does not matter, nor does the order of the keys. The only
mandatory section in any input file is the "[locations]" section. The default
behavior of the `gendef` command is to populate the "[locations]" section with
the path to the current directory. So after creating this file, the only user
work **required** is to finish the paths in the "[locations]" section. The
other sections contain keys that change the default behavior of the pipeline
(**CAUTION**: changing these requires an understanding of the tools involved,
run `lahontan help fcounts` for more much more detailed information)

Additionally, there are certain reserved values that carry particular
significance for the input file parser. The reserved value **None** is reserved
to mean the pipeline will default to its default behavior.  The reserved values
**True** and **False** mean the pipeline will turn on or off the desired
option, respectively.  Other options contain special arguments that you can see
more information about by using the built-in `help` command.

Notice that the default behavior of `gendef` is to default all behavior. Thus,
the previous input file is equivalent to:
```
default_input.ini
-------------------------------------------------------------------------------
[locations]
Project = /home/alberton/rna-seq/scripts/Pipeline/Test
Reference = /home/alberton/rna-seq/scripts/Pipeline/Test
Original = /home/alberton/rna-seq/scripts/Pipeline/Test
```

Finally, a completed input file may look something like:
```
default_input.ini
-------------------------------------------------------------------------------
[locations]
Project = /home/alberton/rna-seq/scripts/Pipeline/Test/MyProject
Reference = /home/alberton/rna-seq/scripts/Pipeline/Test/MyProject-Reference
Original = /home/alberton/rna-seq/scripts/Pipeline/Test/MyProject-Data
```
Notice that the paths to the Project, Reference, and Original directories are
completed.

**Mandatory Arguments**

There are no mandatory arguments for the `gendef` command, however, the `-m`
option specifies the type of input file to create, and it will default to
`fcounts` if not specified.

For example:
```
$ lahontan gendef
```
will create a `fcounts` input file called `default_input.ini`.

Furthermore:
```
$ lahontan gendef -m string -f StringtieInput.ini
```
will create a `string` input file called `StringtieInput.ini`.

#### help
**Usage: lahontan help COMMAND**

Used to display additional information about a command. Similar to `lahontan
<command> --help`

**Available Commands**

```
fcounts     For running featureCounts pipeline
string      For running Stringtie pipeline
kall        For running kallisto pipeline
bowtie2     For running bowtie2 pipeline
mj          Provoke questionnaire to make a Metadata file
mb          Create trimmomatic blacklist
clean       Clean any project directories
fo          Finds optimal execution path for batch execution
```

# Authors

* [**Alberto Nava**](https://github.com/Alberto024)
* [**Richard Tillett**](https://github.com/rltillett)
