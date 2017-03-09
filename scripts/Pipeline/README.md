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
* sleuth
* ballgown
  
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

-e *\<stage\>*, --execute *\<stage\>*  
> Used to specify stage of Pipeline to be ran  
> *\<stage\>* is a comma-separated list of stages to be executed.  
> Possible stages include:  
>> 1 : Creating Project Structure; this includes creating Project directories and making 
symbolic links to reference and raw data  
>> 2 : Preparing Reference Data; this includes running Quality Control then setting up hisat2 databases  
>> 3 : Running actual Pipeline; this includes alignment, compression, and feature counts  
>> 4 : Preparing for R analysis; this includes creating 'nice' count file and developing 
automated R scripts from `Metadata.json`  
>> 5 : Running R analysis; this involves running DESeq2 and edgeR differential expression analyses  
>> A : (1,2,3,4,5) i.e. runs entire pipeline  

> Note: default is to just run `A`

-r *\<sampleNumber\>*, --runsample *\<sampleNumber\>*  
> Used to execute through Stage 3 on the sample specified by *\<sampleNumber\>*  

-s *\<phase>\*, --stringtie *\<phase>\*
> Modifies pipeline to incorporate stringtie and ballgown. Stringtie replaces featureCounts.  
> Possible phases can be any of: "a","b","c","ab","bc","abc"  
>> a : Phase a refers to mapping sample reads to reference genome, converting mapped 
reads to .bam format, and then assembling transcripts for each sample  
>> b : Phase b refers to merging the transcripts from each sample into one merged .gtf file  
>> c : Phase c refers to estimating the transcript abundances for each sample  

> Option also used to specify stringtie postprocessing options in which case use: `--stringtie abc` when executing stages 4 or 5.  
> Note: If you will be running phase b, you cannot specify a sample to run with --runsample

-k, --kallisto  
> Modifies pipeline to incorporate kallisto and sleuth. Kallisto replaces featureCounts and hisat2.  
> Option also used to specify kallisto postprocessing options in which case 
use: `--kallisto` when executing stages 4 or 5.  

-j *\<jsonFile\>*, --jsonfile *\<jsonFile\>*  
> Ignores JSON Metadata file creation and uses specified path to already created JSON Metadata  

-c *\<placeToClean\>*, --clean *\<placeToClean\>*  
> Used to clean/reset *\<placeToClean\>*  
> Possible places include:
>> Reference = removes all files except Genome, cDNA, and gtf in Reference directory  
>> Data = removes all files except raw data file links in Data directory  
>> Postprocessing = removes all files except `Metadata.json` if it exists in Postprocessing directory  
>> All = resets entire roject. This is equivalent to running stage 1 on fresh data

--sampleclean *\<sampleName\>*  
> Used to clean/reset *\<sampleName\>*  
> Similar to `--clean Data` but instead of cleaning all samples just cleans a single sample directory  
> Note: *\<sampleName\>* will most likely be of the form `sample_#` where # refers to specific sample number

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
> Note: If you want to run either the kallisto modified pipeline or the stringtie 
modified pipeline, you must also provide either `--kallisto` or `--stringtie abc` 
somewhere in the command

--makebatchbiox  
> Modifies behavior of `--makebatch`. Makes batch file with best behavior for our cluster (compute-[0-2])  

--batchjson *\<pathtoJSON\>*  
> If an optimal path JSON file has been created using `optPath.py` or another method, 
use option to skip the calculation of the optimal CPU usage  
> Note: If option to be used, must be ran with `--makebatch`  

--noconfirm  
> Used for non-interactivity  
> Ignores all user prompts except JSON file creation unless given `--jsonfile` option  

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
> This will run the pipeline using the input variables from */path/to/INPUT* 
with all the default behavior  

`runPipe --execute 1,2 /path/to/INPUT`
> This will run the first and second stage of the pipeline using the input variables from */path/to/INPUT*  

`runPipe --runsample 1,2 /path/to/INPUT`
> This will run through Stage 3 on samples 1 and 2(sample_01 & sample_02)  

`runPipe --jsonfile /path/to/Metadata /path/to/INPUT`
> This will run the pipeline using the input variables from */path/to/INPUT* and 
will use the Metadata located at */path/to/Metadata*  

`runPipe --clean Data /path/to/INPUT`
> This will wipe all of the pipeline's actions on all samples within the respective Data 
folder in its Project Directory(specified in */path/to/INPUT*) leaving only the 
symbolic links to raw data that were created  
> Note: this leaves the project in a state equivalent to `runPipe --execute 1 /path/to/INPUT`

`runPipe --NUKE /path/to/INPUT`
> This will remove entire Project directory(specified in */path/to/INPUT*)  
> Note: Will leave original Reference and Data folders as they were originally  

```
runPipe --reference-qc /path/to/ReferenceDirectory /path/to/INPUT
runPipe --reference-pp /path/to/ReferenceDirectory /path/to/INPUT
```
> If no errors in reference files are detected by first command, then second 
command can be ran to prepare reference data  
  
`runPipe --makebatch 48,48,32 --batchjson /path/to/Optimization/Slurm/JSON --use-reference --jsonfile /path/to/Metadata --stringtie abc /path/to/INPUT`
> This shows an example of creating a batch file for slurm to be used on a 
> cluster with 3 servers with 48, 48, and 32 CPUs on each. An optimized path has been 
> created using optPath.py and saved to /path/to/Optimization/Slurm/JSON. The reference
> directory has also been prepared using a combination of `--reference-qc`, `--reference-pp`, and
> `--stringtie abc`. A Metadata file has also been created by makeJSON.py and has been
> saved to /path/to/Metadata. Finally, we specify that we wish to use the stringtie modified
> pipeline with `--stringtie abc` and conclude with the INPUT file located at
> /path/to/INPUT.


## Project Structure
  
Stage 1 of the Pipeline will create a consistent directory structure between 
Projects that includes a *Data*, *Original*, *Postprocessing*, and *Reference* directory  
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
  
## Creating a Project and running Pipeline
  
### Requirements  
* Paired-end Illumina fastq files. This means 2 files per sample. One for 
Read1, the other for Read2  
* A reference genome for the organism in study. All chromosomes in a single
 FASTA format file  
* A GTF-format file for containing the gene and exon coordinates for your 
reference genome  
* A reference transcriptome (or some set of known transcripts) for your organism  
* The Pipeline directory should be in your `$PATH` environmental variable. This can 
be accomplished by being in the main rna-seq directory and running:  
```
$ cd scripts/Pipeline
$ echo 'export PATH=$PATH:'`pwd` >> $HOME/.bashrc
$ source $HOME/.bashrc
```
### Setup  
  
Often times your data will be organized in a way that isn't accepted by runPipe. 
In order to continue, your `Original` directory must be arranged such that each 
sample read is chronologically sorted in a single directory (Note: Sorted by 
sample and by read)  

Example Original Directory:  
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

Also, in order to limit errors, your filenames must only contain: letters, 
numbers, periods, hyphens, and underscores.  
Note: `Extract.py` might have some useful functions for organizing files  


Your reference filenames must also follow a few specific rules. Along with only using 
letters, numbers, periods, hyphens, and underscores:  
* Your GTF-format file must end with ".gtf" and it must be the only file within your 
reference directory to do so  
* Your reference transcriptome must contain ".cdna." anywhere within the name except 
for the end. It must also be the only file within your reference directory to do so  
* Your reference genome must contain ".dna." anywhere within the name except for the end.
 It must also be the only file within your reference directory to do so  

Example Reference Directory:  
```
Blah_Plantia-Reference/
    Blah_Plantia.release-10.cdna.all.fa.gz
    Blah_Plantia.release-10.dna.toplevel.fa.gz
    Blah_Plantia.release-10.gtf
```

#### Optional but Recommended Setup

I recommend you create a `Metadata.json` file before starting the pipeline if you want 
R analyses. The advantage of creating it beforehand is that you won't be spending 
valuable CPU time and you'll be able to edit it if you mess up during the questionnaire. 
To create your Metadata file, there is a script `makeJSON.py` which asks you a series of 
questions asking you to identify features of your experiment.  
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

I also recommend you prepare your reference data before running the pipeline simply 
because it doesn't require as many resources as some of the other tools require, so 
it would be ideal to prepare it beforehand which also gives you the privilege of time 
as you won't be wasting valuable CPU resources.  

In order to prepare your reference data, you will most likely need to have some knowledge of fasta format.  
To begin, in the best case that your reference data is perfect, just run:
```
$ runPipe --reference-qc /path/to/ReferenceDirectory /path/to/INPUT
$ runPipe --reference-pp /path/to/ReferenceDirectory /path/to/INPUT
```
Note: the INPUT file is a required argument but since it isn't being read, you can use any ordinary file  

In some cases your reference genome will contain revised chromosomes, but since we want 
only the true chromosomes, you can use:  
```
$ seqtk subseq Blah_Plantia.release-10.dna.toplevel.fa.gz Chromosomes_We_Want.txt > Blah_Plantia.TrueChromosomes.dna.fa.gz
```

In other cases the reference transcriptome doesn't agree with the GTF-format file. 
In that case, you can generate a new cDNA file by using:  
```
$ gffread Blah_Plantia.release-10.gtf -g Blah_Plantia.TrueChromosomes.dna.fa.gz -w Blah_Plantia.ours.cdna.fa
```

Lastly, if you will be using slurm for your computation I recommend you create an optimal 
path file using `optPath.py`. The optimal path file is used to determine the best way to 
allocate resources for each of the jobs within the pipeline. The algorithm that I use is 
very robust but under certain paramaters has poor performance. For this reason, I recommend 
you create an optimal path file before starting the pipeline.  

`optPath.py` takes 3 command line arguments. These are:  
1. The number of samples in your experiment  
2. The number of processors to be used to calculate optimal path  
3. A comma-separated list of CPUs on each node in your cluster  
For example:  
```
$ optPath.py 4 48 48,48,32
```
This is telling `optPath.py` to calculate using 48CPUs the best way to run 4 samples 
on a cluster containing two 48CPU computers as well as a 32CPU computer. `optPath.py` 
will create a file called `OptimalPath.dat`:  
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
This says that on a cluster of three nodes that has 48CPU on two nodes and a third with 
32CPU, the optimal usage is to run one sample on each 48CPU node while leaving the 32CPU 
node unused and then doing the same thing for the other two samples. The pipeline will 
use this information to tell slurm how to optimally distribute the jobs.  
Note:  `optPath.py` might require some tinkering if convergence times are slow.  

## Running the Pipeline with slurm

If you will be using slurm to run your pipeline the next step is to create your batch 
file. If you have followed the recommended setup then the command used to create 
your batch file should be:  

```
$ runPipe --jsonfile /path/to/Metadata.json --use-reference --makebatch 48,48 --batchjson /path/to/OptimalPath.dat /path/to/INPUT
```

Note that the argument to `--makebatch` was "48,48" instead of "48,48,32". Because we found out from our 
`OptimalPath.dat` that only the two biggest nodes would be used so we can leave the third smaller node 
unallocated. We also used the Metadata file that we created by passing `--jsonfile` option. We passed 
the `--batchjson` option to use the `OptimalPath.dat` file that tells us the best way to allocate the 
jobs. Make sure to pass the `--use-reference` argument if you prepared your reference data beforehand. 

Upon successful creation of your batch file, you should recieve a message that looks like:
```
$ runPipe --jsonfile /path/to/Metadata.json --use-reference --makebatch 48,48 --batchjson /path/to/OptimalPath.dat /path/to/INPUT
Batch file successfully created:
        /path/to/pipeBatch
```

This will give you a batch file that should look something like this:
```
/path/to/pipeBatch
---------------------------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=400
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=2
#SBATCH --job-name="Pipeline"
#SBATCH --export=PATH,RNASEQDIR,HOME

inputFile='/path/to/INPUT'
jsonFile='/path/to/Metadata.json'

# Stage 1 and 2
srun -N1 -c1 -n1 runPipe --noconfirm --use-reference --jsonfile "${jsonFile}" --execute 1,2 "${inputFile}"

wait

# Stage 3
srun -N1 -c48 -n1 --exclusive runPipe --noconfirm --use-reference --maxcpu 48 -e 3 -r 1 "${inputFile}" &
srun -N1 -c48 -n1 --exclusive runPipe --noconfirm --use-reference --maxcpu 48 -e 3 -r 2 "${inputFile}" &
wait
srun -N1 -c48 -n1 --exclusive runPipe --noconfirm --use-reference --maxcpu 48 -e 3 -r 3 "${inputFile}" &
srun -N1 -c48 -n1 --exclusive runPipe --noconfirm --use-reference --maxcpu 48 -e 3 -r 4 "${inputFile}" &

wait

# Stage 4
srun -N1 -c1 -n1 runPipe --noconfirm --use-reference --jsonfile "${jsonFile}" --execute 4 "${inputFile}"

wait

# Stage 5
srun -N1 -c1 -n1 --exclusive runPipe --noconfirm --use-reference --jsonfile "${jsonFile}" --execute 5 --edger "${inputFile}" &
srun -N1 -c1 -n1 --exclusive runPipe --noconfirm --use-reference --jsonfile "${jsonFile}" --execute 5 --deseq "${inputFile}" &

scontrol show job $SLURM_JOB_ID
wait
```
You should probably look through it to look for any errors and edit any parameters. 
For example, you might wish to edit the `--job-name` that slurm will give your job. 
Also, if you do not wish to run the R analyses that I have provided, then you can 
comment out both the Stage 5 commands leaving only the last 2 lines at the bottom of 
the file. Once you are satisfied with your batch, you can run it with:
```
$ sbatch /path/to/pipeBatch &
```
  
### Running a Stringtie Modified Pipeline with slurm

If you wish to run a stringtie modified pipeline with slurm, the only distinction that
you need to make is to add `--stringtie abc` to your command that creates your batch
file:  
  
```
$ runPipe --jsonfile /path/to/Metadata.json --use-reference --makebatch 48,48 --batchjson /path/to/OptimalPath.dat --stringtie abc /path/to/INPUT
```
  
The main difference within the created batch file should be that your commands have a
corresponding `--stringtie` options, but you should still throroughly look over it.

### Running a Kallisto Modified Pipeline with slurm

If you wish to run a kallisto modified pipeline with slurm, the only distinction that
you need to make is to add `--kallisto` to your command that creates your batch
file:  
  
```
$ runPipe --jsonfile /path/to/Metadata.json --use-reference --makebatch 48,48 --batchjson /path/to/OptimalPath.dat --kallisto /path/to/INPUT
```
  
The main difference within the created batch file should be that your commands have a
corresponding `--kallisto` options, but you should still thoroughly look over it.
  

## Running the Pipeline without slurm

If you don't have access to a cluster fitted with slurm or you just wish to test the 
pipeline on a more interactive scale then `runPipe` can be ran just as well.  

To begin, if you have followed the proper setup, then you should have:  
* An `Original` directory containing your raw data files sorted in the proper order  
* A `Reference` directory containing your genome, reference transcriptome, and GTF 
file along with the other files that were prepared using `runPipe --reference-pp`  
* An INPUT file that defines paths to those `Original` and `Reference` directories 
as well as a `Project` path  
* A correct and valid Metadata file in JSON format created by `makeJSON.py`  

Once you have those four things you can create your Project structure by running:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 1,2 /path/to/INPUT
```
Note that the default behavior of `runPipe` is to complete its tasks as quickly as 
possible. If you wish to leave some CPU resources unallocated, which I would recommend 
you do, you can pass `--maxcpu` to `runPipe`  

You can find out how many CPU your computer has by running:  
```
$ nproc
```
So if `nproc` returns 4 and we want to leave 1 CPU open, the previous command becomes:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 1,2 --maxcpu 3 /path/to/INPUT
```
Once you have created your Project structure you can run the most CPU-intensive Stage 3:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 3 --maxcpu 3 /path/to/INPUT
```
The default behavior is to run every sample but if you wish to run any give sample, 
just pass `--runsample`. For example, if you only wish to run the first sample, 
the previous command becomes:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 3 --maxcpu 3 --runsample 1 /path/to/INPUT
```
Once every sample has completed, you can run Stage 4 which prepares R scripts as 
well as creating a "Nice" feature counts file:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 4 --maxcpu 3 /path/to/INPUT
```
There will be a `Counts.dat`,`NiceCounts.dat`, and a `GoodCounts.dat` located in 
the `Postprocessing` directory of your Project.
* `Counts.dat` will give you the least information. It has only the ID of the gene 
and its respective counts  
* `NiceCounts.dat` has the ID of the gene, the length of the gene, as well as 
its respective counts  
* `GoodCounts.dat` will give you the most information. It has the gene ID, gene Name, 
gene biotype, chromosome location, length, and its respective counts.  

If you do not wish to run R analyses, the hopefully the Pipeline went smoothly 
and were happy with the results. If you do wish to run R analyses you can run:
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 5 --maxcpu 3 /path/to/INPUT
```
or
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 5 --maxcpu 3 --deseq /path/to/INPUT
```
or
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 5 --maxcpu 3 --edger /path/to/INPUT
```
The first command will run both DESeq2 and edgeR analyses. The second command will run 
only DESeq2 and the third will run only edgeR. The DESeq2 analysis will create simple 
and more detailed PCA plots as well as creating reports using `regionReport` and 
`ReportingTools`. You can view what the DESeq2 analysis is doing by looking at 
`makeReport.r`. The edgeR analysis will create a report using `regionReport`. You can 
view what the edgeR analysis is doing by looking at `makeEdge.r`.

### Running Stringtie Modified Pipeline without slurm

Running the Stringtie Modified Pipeline is very similar to running the typical
featureCounts pipeline, however, there are a few key differences. Stages 1 and 2 should
be ran with the same command with an additional `--stringtie abc` option. This
is to inform the pipeline to process the data with Stringtie. So stages 1 and 2 should 
look like:  
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 1,2 --stringtie abc /path/to/INPUT
```
Stage 3 is where the biggest differences occur. Stage 3 for the stringtie modified
pipeline is split up into 3 phases: phase a, phase b, and phase c.   
* Phase a refers to mapping each of the sample reads to a reference genome, converting the mapped reads
to bam format with samtools, and then assembling the transcripts for each sample. These
are all steps that occur within the featureCounts pipeline as well. Up through phase a
the featureCounts and stringtie pipelines are more or less identical.   
* Phase b of the third
stage within the stringtie pipeline is where the analysis diverges. Phase b refers to
merging the transcripts from each sample into one merged gene transcript file.  
* Phase c refers to estimating the transcript abundances for each sample.  

Phase a and c are processes that
run on each sample individually, while phase b runs on the entire data set -- 
meaning **phase a and c can be ran with a `--runsample` option while phase b cannot**.

The argument to the `--stringtie` option corresponds to the phase of stage 3 that you
want to run, except when denoting any other stage, then you must use `--stringtie abc`.  

So, if you want to run stage 3 in its entirety, you should use `--stringtie abc` which will run
phases a, b, and c in sequence. If you just want to run phase a, then use `--stringtie a`.
If you want to run phase b, then use either `--stringtie ab` or `--stringtie b` but note that 
phase b will not be able to run if phase a has not been completed, so it is recommended to run
the former option.  

Finally, here is an example:  
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 3 --maxcpu 3 --stringtie abc /path/to/INPUT
```

Stage 4 is more or less similar to the normal pipeline with the exception that you must add
the `--stringtie abc` option.  
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 4 --maxcpu 3 --stringtie abc /path/to/INPUT
```

Similarly for Stage 5, but note that stringtie has only one R analysis option, ballgown,
so there are no extra available options for Stage 5 besides `--stringtie abc`.  
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 5 --maxcpu 3 --stringtie abc /path/to/INPUT
```

### Running Kallisto Modified Pipeline without slurm

Running the Kallisto Modified Pipeline is also very straightforward. Unlike
the Stringtie pipeline, the Kallisto pipeline has only one phase during stage 3
so it can be run identically as the normal pipeline with the exception that
you must add the `--kallisto` option for each command. So the entire pipeline
may be ran with:   
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 1,2 --kallisto /path/to/INPUT
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 3 --maxcpu 3 --kallisto /path/to/INPUT
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 4 --maxcpu 3 --kallisto /path/to/INPUT
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --execute 5 --maxcpu 3 --kallisto /path/to/INPUT
```
or in one command:   
```
$ runPipe --use-reference --jsonfile /path/to/Metadata.json --kallisto /path/to/INPUT
```
