'Usage: runSleuth.r [-h | --help] [-p </path/to/data>] [--cols <name>] [--base <name>]

Example:
Rscript runSleuth.r -p /path/to/project --cols ColumnsName

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols file [default: ./]
    --cols <name>       Name of phenotype file that describes experiment [default: Cols.dat]
    --base <name>       Name of Base Directory that contains Kallisto ouput 
                        [default: KallistoResults]'-> doc

# Set up command line argument library
library('docopt')

# Retrieve command line arguments
opts <- docopt(doc)

# Change to directory
setwd(opts$p)

# Load the necessary libraries
library('devtools')
library('sleuth')

base_dir <- opts$p
sample_id <- dir(file.path(base_dir,opts$base))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, opts$base, id, "abundance.h5"))
dataFrame <- read.csv(opts$cols, header=TRUE, stringsAsFactors=FALSE)
dataFrame <- dplyr::mutate(dataFrame, path=kal_dirs)
dataFrame
model <- sleuth_prep(dataFrame, ~ condition)
model <- sleuth_fit(model)
model <- sleuth_fit(model, ~1, 'reduced')
model <- sleuth_lrt(model, 'reduced', 'full')
models(model)
#sleuth_live(model)
results_table <- sleuth_results(model, 'reduced:full', test_type='lrt')
write.csv(results_table, "sleuthKallistoAnalysis.csv", row.names=FALSE)

# unload the libraries
detach('package:docopt')
detach('package:sleuth')
detach('package:devtools')
