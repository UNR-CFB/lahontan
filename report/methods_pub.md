# Methods & Materials
## Library preparation
[Screen 1]
[Screen 2]

## Illumina sequencing
Transcript libraries were sequenced using multiplexed Illumina HiSeq 2000 sequencing, performed by sequencing provider BGI.
[Details of sequencing as performed at BGI should go here. Did they provide boiler-plate?]

 Paired-end sequencing was performed at 2 x 50 bp length for the [screen1] samples. [Screen2] libraries were sequenced at 2 x 100 bp.

## Sequence quality control
Sequence pairs were trimmed and filtered for nucleotide-base quality and sequencing adapters using Trimmomatic v. 0.35 (Bolger et al, 2014), with settings as follows: removal of Illumina TruSeq paired-end adapter sequences, trimming of low-quality  (phred quality ≤ 5) and/or “N” bases from both ends of each read, further trimming using a 4-base sliding window to remove regions with average Q < 5, and, finally, removal of trimmed sequences with length < 25 nt. Sequence quality was visualized before and after trimming using FastQC v. 0.11.4 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).


## Sequence alignment and expression quantification
Trimmed sequence pairs were aligned to the Ensembl reference mouse genome (Ensembl build 83; GRCm38.p4) using the HISAT2 spliced read alignment tool (v. 2.0.1; Kim et al, 2015) [Note: pub is for HISAT1, HISAT2 not yet pub’d]. For Mus musculus genes annotated in Ensembl build 83 (Yates et al, 2016) in Gene Transfer Format (GTF), the genomic coordinates of all transcriptional splice sites were extracted and supplied to HISAT2 at time of alignment (via the `--known-splicesite-infile` option). Aligned outputs were converted to the binary BAM format (Li et al, 2009).

Following alignment, the raw counts of read pairs aligned to each gene were totalled using the featureCounts tool of the subread package (v. 1.5.1-p1; Liao et al, 2014). Reads were counted just once per pair, summarized for gene loci, with only read pairs aligned to a single transcribed location included in the count totals.

## Quantitative Analysis
[KAREN SECTION]
Count data underwent low-count filtering, upper quartile normalization, and standard quality control before they were transformed into RPKM.



# References
Bolger AM, Lohse M, Usadel B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, 30(15):2114-2120.
Kim D, Langmead B, Salzberg SL (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360.
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. (2009) The sequence alignment/map format and SAMtools. Bioinformatics. 25(16):2078-9.
Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30.
Yates A, Akanni W, Amode MR, Barrell D, Billis K, Carvalho-Silva D, Cummins C, Clapham P, Fitzgerald S, Gil L, Girón CG, Gordon L, Hourlier T, Hunt SE, Janacek SH, Johnson N, Juettemann T, Keenan S, Lavidas I, Martin FJ, Maurel T, McLaren W, Murphy DN, Nag R, Nuhn M, Parker A, Patricio M, Pignatelli M, Rahtz M, Riat HS, Sheppard D, Taylor K, Thormann A, Vullo A, Wilder SP, Zadissa A, Birney E, Harrow J, Muffato M, Perry E, Ruffier M, Spudich G, Trevanion SJ, Cunningham F, Aken BL, Zerbino DR, Flicek P. (2016) Ensembl 2016. Nucleic acids research, 44(D1):D710-6.
