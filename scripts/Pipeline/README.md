# runPipe
- - -
  
## Dependencies
* python3  
* python3-docopt
* bash
* shutils
  
### Optional (used for R analysis)
* R >= 3.3
* DESeq2
* edgeR
* docopt
* ggplot2
* ReportingTools
* regionReport
  
## Usage: runPipe [options] *\<pathToInputFile\>*
  
## Arguments
*\<pathToInputFile\>* = a valid path to an input file  
The input file should contain definitions for three paramaters:  
* Project = Directory where Pipeline should be ran and saved to  
* Reference = Directory containing reference materials  
** Reference directory should contain at least: a Genome, cDNA, and .gtf file
* Original = Directory containing raw data files  

Example Input file:
```
sampleInputFile
---------------------------------------------------------------------------------------------------
Project="/home/user/PlantExperiment"
Reference="/home/user/Reference/Plant"
Original="/home/user/Data/PlantTestA"
```
  
## Options
  
-h, --help  
> Shows usage and options, then exits

--version  
> Shows `runPipe` version number, then exits

-j *\<jsonFile\>*, --jsonfile *\<jsonFile\>*  
> Ignores JSON Metadata file creation and uses specified path to already created JSON Metadata  

-c *\<placeToClean\>*, --clean *\<placeToClean\>*  
> Used to clean/reset *\<placeToClean\>*  
> Possible places include:
>> Reference = removes all files except Genome, cDNA, and gtf in Reference directory  
>> Data = removes all files except raw data file links in Data directory  
>> Postprocessing = removes all files except `Metadata.json` if it exists in Postprocessing directory  
>> All = resets entire roject. This is equivalent to running stage 1 on fresh data

-s *\<sampleName\>*, --sampleclean *\<sampleName\>*  
> Used to clean/reset *\<sampleName\>*  
> Similar to `--clean Data` but instead of cleaning all samples just cleans a single sample directory  
> Note: *\<sampleName\>* will most likely be of the form `sample_#` where # refers to specific sample number

-e *\<stage\>*, --execute *\<stage\>*  
> Used to specify stage of Pipeline to be ran  
> *\<stage\>* is a comma-separated list of stages to be executed.  
> Possible stages include:  
>> 1 : Creating Project Structure; this includes creating Project directories and making symbolic links to reference and raw data  
>> 2 : Preparing Reference Data; this includes running Quality Control then setting up hisat2 databases  
>> 3 : Running actual Pipeline; this includes alignment, compression, and feature counts  
>> 4 : Preparing for R analysis; this includes creating 'nice' count file and developing automated R scripts from `Metadata.json`  
>> 5 : Running R analysis; this involves running DESeq2 and edgeR differential expression analyses  
>> A : (1,2,3,4,5) i.e. runs entire pipeline  

> Note: default is to just run `A`

-r *\<sampleNumber\>*, --runsample *\<sampleNumber\>*  
> Used to execute through Stage 3 on the sample specified by *\<sampleNumber\>*  

--maxcpu *\<CPUs\>*  
> Used to limit number of CPUs used by Pipeline  
> Note: default is to use all available CPUs  

--reference-qc *\<pathtoReferenceDirectory\>*  
> Used to prepare reference data outside of the project structure  
> Runs Quality Control check on Reference files i.e. First part of Stage 2  

--reference-pp *\<pathtoReferenceDirectory\>*  
> Used to prepare reference data outside of the project structure  
> Pre-processes Reference data i.e. Second part of Stage 2  

--use-reference  
> If reference data has already been prepared, use option to skip Stage 2 and use prepared reference data  
> Note: If option is to be used, must be ran from Stage 1  

--makebatch *\<cluster\>*  
> Used to make batch file for use with slurm  
> *\<cluster\>* is a comma-separated list of CPUs on each node in your cluster  

--makebatchbiox  
> Modifies behavior of `--makebatch`. Makes batch file with best behavior for our cluster (compute-[0-2])  

--batchjson *\<pathtoJSON\>*  
> If an optimal path JSON file has been created using `optPath.py` or another method, use option to skip the calculation of the optimal CPU usage  
> Note: If option to be used, must be ran with `--makebatch`  

--noconfirm  
> Used for non-interactivity  
> Ignores all user prompts except JSON file creation  

--NUKE  
> Removes entire project  

--edger  
> Runs edgeR analysis only  
> Note: default is to run both  

--deseq  
> Runs DESeq2 analysis only  
> Note: default is to run both  

## Examples
  
`runPipe --help`
> Will show you available options

`runPipe /path/to/INPUT`
> This will run the pipeline using the input variables from */path/to/INPUT* with all the default behavior  

`runPipe --execute 1,2 /path/to/INPUT`
> This will run the first and second stage of the pipeline using the input variables from */path/to/INPUT*  

`runPipe --runsample 1,2 /path/to/INPUT`
> This will run through Stage 3 on samples 1 and 2(sample_01 & sample_02)  

`runPipe --jsonfile /path/to/Metadata /path/to/INPUT`
> This will run the pipeline using the input variables from */path/to/INPUT* and will use the Metadata located at */path/to/Metadata*  

`runPipe --clean Data /path/to/INPUT`
> This will wipe all of the pipeline's actions on all samples within the respective Data folder in its Project Directory(specified in */path/to/INPUT*) leaving only the symbolic links to raw data that were created  
> Note: this leaves the project in a state equivalent to `runPipe --execute 1 /path/to/INPUT`

`runPipe --NUKE /path/to/INPUT`
> This will remove entire Project directory(specified in */path/to/INPUT*)  
> Note: Will leave original Reference and Data folders as they were originally  

```
runPipe --reference-qc /path/to/ReferenceDirectory /path/to/INPUT
runPipe --reference-pp /path/to/ReferenceDirectory /path/to/INPUT
```
> If no errors in reference files are detected by first command, then second command can be ran to prepare reference data  
  
## Project Structure
  
Stage 1 of the Pipeline will create a consistent directory structure between Projects that includes a *Data*, *Original*, *Postprocessing*, and *Reference* directory  
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
* *Data* = directory where Stage 3 will propogate. Contains a directory for each respective sample which contains its respective raw data files  
* *Original* = directory containing symbolic links pointing to original raw data files  
* *Postprocessing* = directory where differential expression analysis and other analyses will be saved to  
* *Reference* = directory containing symbolic links pointing to original reference files  
  
## Creating a Project and running Pipeline
  
### Requirements  
* Paired-end Illumina fastq files. This means 2 files per sample. One for Read1, the other for Read2. The examples will use files named read1.fq.gz and read2.fq.gz  
* A reference genome for the organism in study. All chromosomes in a single FASTA format file  
* A GTF-format file for containing the gene and exon coordinates for your reference genome  
* A reference transcriptome (or some set of known transcripts) for your organism  
* The Pipeline directory should be in your `$PATH` environmental variable. This can be accomplished by being in the main rna-seq directory and running:  
```
cd scripts/Pipeline
echo 'export PATH=$PATH:'`pwd` >> $HOME/.bashrc
source $HOME/.bashrc
```
### Setup  
  
Often times your data will be organized in a way that isn't accepted by runPipe. In order to continue, your raw data must be arranged such that each sample read is chronologically sorted in a single directory (Note: Sorted by sample and by read)  

Example Raw Data Directory:  
```
Blah_Plantia-RawData/
    Blah_Plantia.treament1_group1.read1.fq.gz
    Blah_Plantia.treament1_group1.read2.fq.gz
    Blah_Plantia.treament1_group2.read1.fq.gz
    Blah_Plantia.treament1_group2.read2.fq.gz
    Blah_Plantia.treament2_group1.read1.fq.gz
    Blah_Plantia.treament2_group1.read2.fq.gz
    Blah_Plantia.treament2_group2.read1.fq.gz
    Blah_Plantia.treament2_group2.read2.fq.gz
```

Also, in order to limit errors, your filenames must only contain: letters, numbers, periods, hyphens, and underscores.  
Note: `Extract.py` might have some useful functions for organizing files  


Your reference filenames must also follow a few specific rules. Along with only using letters, numbers, periods, hyphens, and underscores:  
* Your GTF-format file must end with ".gtf" and it must be the only file within your reference directory to do so  
* Your reference transcriptome must contain ".cdna." anywhere within the name except for the end. It must also be the only file within your reference directory to do so  
* Your reference genome must contain ".dna." anywhere within the name except for the end. It must also be the only file within your reference directory to do so  

Example Reference Directory:  
```
Blah_Plantia-Reference/
    Blah_Plantia.release-10.cdna.all.fa.gz
    Blah_Plantia.release-10.dna.toplevel.fa.gz
    Blah_Plantia.release-10.gtf
```

#### Optional but Recommended Setup

I recommend you create a `Metadata.json` file before starting the pipeline if you want R analyses. The advantage of creating it beforehand is that you won't be spending valuable CPU time and you'll be able to edit it if you mess up during the questionnaire. To create your Metadata file, there is a script `makeJSON.py` which asks you a series of questions asking you to identify features of your experiment.  
Note: Your features should contain at least one letter.  
```
$ makeJSON.py
What is the name of the project? Blah_Plantia-Experiment  
How many samples are there? 4
How many features does each sample have? 2
What is the name of feature #1? treatment
What is the name of feature #2? group
What is the name of the main feature? treatment
What is the treatment feature for sample #1? one
What is the group feature for sample #1? one
What is the treatment feature for sample #2? one
What is the group feature for sample #2? two
What is the treatment feature for sample #3? two
What is the group feature for sample #3? one
What is the treatment feature for sample #4? two
What is the group feature for sample #4? two
Done making Metadata.json
```
This creates:
```
Metadata.json
---------------------------------------------------------------------------------------------------
{
    "FeatureNames": [
        "group",
        "treatment"
    ],
    "MainFeature": "treatment",
    "NumberofFeatures": 2,
    "NumberofSamples": 4,
    "ProjectName": "Blah_Plantia-Experiment",
    "Samples": {
        "sample_01": {
            "Features": {
                "group": "one",
                "treatment": "one"
            }
        },
        "sample_02": {
            "Features": {
                "group": "two",
                "treatment": "one"
            }
        },
        "sample_03": {
            "Features": {
                "group": "one",
                "treatment": "two"
            }
        },
        "sample_04": {
            "Features": {
                "group": "two",
                "treatment": "two"
            }
        }
    }
}
```
... which will later be scraped in Postprocessing.  

I also recommend you prepare your reference data before running the pipeline simply because it doesn't require as many resources as some of the other tools require, so it would be ideal to prepare it beforehand which also gives you the privilege of time as you won't be wasting valuable CPU resources.  

In order to prepare your reference data, you will most likely need to have some knowledge of fasta format.  
To begin, in the best case that your reference data is perfect, just run:
```
runPipe --reference-qc /path/to/ReferenceDirectory /path/to/INPUT
runPipe --reference-pp /path/to/ReferenceDirectory /path/to/INPUT
```
Note: the INPUT file is a required argument but since it isn't being read, you can use any ordinary file  

In some cases your reference genome will contain revised chromosomes, but since we want only the true chromosomes, you can use:  
```
seqtk subseq Blah_Plantia.release-10.dna.toplevel.fa.gz Chromosomes_We_Want.txt > Blah_Plantia.TrueChromosomes.dna.fa.gz
```

In other cases the reference transcriptome doesn't agree with the GTF-format file. In that case, you can generate a new cDNA file by using:  
```
gffread Blah_Plantia.release-10.gtf -g Blah_Plantia.TrueChromosomes.dna.fa.gz -w Blah_Plantia.ours.cdna.fa
```

Lastly, if you will be using slurm for your computation I recommend you create an optimal path file using `optPath.py`. The optimal path file is used to determine the best way to allocate resources for each of the jobs within the pipeline. The algorithm that I use is very robust but under certain paramaters has poor performance. For this reason, I recommend you create an optimal path file before starting the pipeline.  

`optPath.py` takes 3 command line arguments. These are:  
1. The number of samples in your experiment  
2. The number of processors to be used to calculate optimal path  
3. A comma-separated list of CPUs on each node in your cluster  
For example:  
```
optPath.py 4 48 48,48,32
```
This is telling `optPath.py` to calculate using 48CPUs the best way to run 4 samples on a cluster containing two 48CPU computers as well as a 32CPU computer. `optPath.py` will create a file called `OptimalPath.dat`:  
```
OptimalPath.dat
---------------------------------------------------------------------------------------------------
{
    "Step 1": {
        "Procs": 48,
        "Samps": 2
    },
    "Step 2": {
        "Procs": 48,
        "Samps": 2
    }
}
```
This says that on a cluster of three nodes that has 48CPU on two nodes and a third with 32CPU, the optimal usage is to run one sample on each 48CPU node while leaving the 32CPU node unused and then doing the same thing for the other two samples. The pipeline will use this information to tell slurm how to optimally distribute the jobs.  
Note:  `optPath.py` might require some tinkering if convergence times are slow.  

## Running the Pipeline with slurm

## Running the Pipeline without slurm

