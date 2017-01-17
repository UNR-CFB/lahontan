# runPipe
***

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

## Usage
### runPipe [options] \<pathToInputFile\>
Note: the Pipeline directory should be in your `$PATH` environmental variable
This can be accomplished by being in the main rna-seq directory and running:
```
cd scripts/Pipeline
echo 'export PATH=$PATH:'`pwd` >> $HOME/.bashrc
source $HOME/.bashrc
```

## Arguments
\<pathToInputFile\> = a valid path to an input file  
The input file should contain definitions for three paramaters:
  * Project = Location of where Pipeline should be ran and saved to
  * Reference = Location of reference materials
  * Original = Location of raw data files  

Example Input file:
```
sampleInputFile
------------------------------------------------------------------------------------------------------------
Project="/home/user/PlantExperiment"
Reference="/home/user/Reference/Plant"
Original="/home/user/Data/PlantTestA"
```
  
## Options
