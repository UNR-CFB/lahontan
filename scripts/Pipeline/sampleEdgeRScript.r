'Usage: makeEdge.r [-h | --help] [-p </path/to/data>] [--counts <name>] [--cols <name>]

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols.dat and Counts.dat [default: ./]
    --counts <name>     Name of collected featureCounts data file [default: Counts.dat]
    --cols <name>       Name of Column file that describes Counts [default: Cols.dat]

Example:
Rscript makeEdge.r -p /path/to/project --counts CountName --cols ColumnsName

(Or if in directory with "Counts.dat" and "Cols.dat", can just call with:
    Rscript makeEdge.r
)'-> doc

# Set up command line argument library
library('docopt')

# Retrieve Command Line Arguments
opts <- docopt(doc)

# change to the directory
setwd(opts$p)

# load the necessary libraries
library('edgeR')
library('ggplot2')
library('regionReport')

# Load the datasets
countData <- read.csv(opts$counts,sep='\t',header=TRUE)
colData <- read.csv(opts$cols,sep='\t',header=TRUE)

# Make initial data frame
y <- DGEList(counts=countData)

# Format colData for edgeR
Group <- factor(paste(colData$condition, sep="."))
cbind(colData, Group=Group)
design <- model.matrix(~ 0+Group, data=y$samples)
colnames(design) <- levels(Group)

# Make edgeR Data Frame
y <- DGEList(counts=countData, group=Group)
y <- calcNormFactors(y) 

# Filter out rows that contain only 1 or 0
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# Calculating Dispersion
y <- estimateDisp(y, design)

################################################################
# Actual Analysis
################################################################

###########################
# Quasi-likelihood F-test #
###########################
# Fit for F-test
QLfit <- glmQLFit(y, design)
# Making contrast matrix for QL F-test
QLcon <- makeContrasts(
	Wtctrl_vs_cngcHS = Wtctrl - cngcHS,
	Wtctrl_vs_WtHS = Wtctrl - WtHS,
	Wtctrl_vs_cngcctrl = Wtctrl - cngcctrl,
	cngcHS_vs_WtHS = cngcHS - WtHS,
	cngcHS_vs_cngcctrl = cngcHS - cngcctrl,
	WtHS_vs_cngcctrl = WtHS - cngcctrl, levels=design)
# QL F-test
anova <- glmQLFTest(QLfit, contrast=QLcon)
topTags(anova)

#############
# GLM Stuff #
#############
# Generalized Linear Model Fit
GLfit <- glmFit(y, design)

# Perform all possible comparisons
mycontrasts <- makeContrasts(
	WtctrlvsWtHS = Wtctrl - WtHS,
	WtctrlvscngcHS = Wtctrl - cngcHS,
	Wtctrlvscngcctrl = Wtctrl - cngcctrl,
	cngcHSvsWtHS = cngcHS - WtHS,
	WtHSvscngcctrl = WtHS - cngcctrl,
	cngcHSvscngcctrl = cngcHS - cngcctrl,
	levels = design)

lrt_WtctrlvsWtHS <- glmLRT(GLfit, contrast=mycontrasts[,"WtctrlvsWtHS"])
lrt_WtctrlvscngcHS <- glmLRT(GLfit, contrast=mycontrasts[,"WtctrlvscngcHS"])
lrt_Wtctrlvscngcctrl <- glmLRT(GLfit, contrast=mycontrasts[,"Wtctrlvscngcctrl"])
lrt_cngcHSvsWtHS <- glmLRT(GLfit, contrast=mycontrasts[,"cngcHSvsWtHS"])
lrt_WtHSvscngcctrl <- glmLRT(GLfit, contrast=mycontrasts[,"WtHSvscngcctrl"])
lrt_cngcHSvscngcctrl <- glmLRT(GLfit, contrast=mycontrasts[,"cngcHSvscngcctrl"])

topTags(lrt_WtctrlvsWtHS)
topTags(lrt_WtctrlvscngcHS)
topTags(lrt_Wtctrlvscngcctrl)
topTags(lrt_cngcHSvsWtHS)
topTags(lrt_WtHSvscngcctrl)
topTags(lrt_cngcHSvscngcctrl)

# Make regionReport
report <- edgeReport(dge = y, object = lrt_WtctrlvsWtHS,
                       project='edgeR Analysis',
                       intgroup='group',
                       outdir='edgeR_Report',
                       ouput='index',
                       theme=theme_linedraw(),
                       output_format='pdf_document',
                       device='pdf')

# unload the libraries
detach("package:edgeR")
detach("package:ggplot2")
detach("package:regionReport")
detach("package:docopt")
