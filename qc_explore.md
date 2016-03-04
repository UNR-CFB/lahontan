# QC and alignment stats walkthrough
## Big picture
_1-2 sentence version:_ screen 1 (8 samples, 4x2) looked absolutely great, from qc to counts; screen 2 (2 x 2, low-input RNA) samples contain  artificial/artifact sequences we believe originate from their more challenging preps, which drop our alignment success rates, but otherwise should not harm results.

Below, we'll dig into what those steps involved and found, starting with screen 1:
## Screen 1 QC and alignment stats
### QC (fastQC)
For screen1, the 8 sample (4x2) RNA-seq experiment, the QC looked typical and good for Illumina RNA-seq.
QC data (Fastqc) looked good, typical for RNA-seq and overall a-OK.
Fastqc reports about a dozen results in a handy html file. Example: here's the output for from the read1 (of 2, paired-end) sequences of sample1:
* http://bioinformatics.unr.edu/share/jr/gould_qc/read_1_fastqc.html
* Basic info: ~41M reads, 100 letters long per alf, Illumina 1.5 (phred64) file encoding
* The per-base sequence content warning is common with Illumina reads, and likely relates to biases in the Illumina chemistry.
* The sequence duplication levels are also frequently seen in Illumina RNA-seq data.
* Overall everything looks a-OK.

### Alignment (HISAT2)
After trimming and re-running fastqc on the trimmed reads, we aligned all 41M pairs to the mouse genome (and all known splice sites) with HISAT2, and the result summary looked great:
```
41148480 reads; of these:
  41148480 (100.00%) were paired; of these:
    4407760 (10.71%) aligned concordantly 0 times
    30601389 (74.37%) aligned concordantly exactly 1 time
    6139331 (14.92%) aligned concordantly >1 times
    ----
    4407760 pairs aligned concordantly 0 times; of these:
      394864 (8.96%) aligned discordantly 1 time
    ----
    4012896 pairs aligned 0 times concordantly or discordantly; of these:
      8025792 mates make up the pairs; of these:
        5828721 (72.62%) aligned 0 times
        1602130 (19.96%) aligned exactly 1 time
        594941 (7.41%) aligned >1 times
92.92% overall alignment rate
```
* 30M pairs both aligned exactly once -- fantastic.
* 2 M solo reads aligned exactly once (partner read aligned zero times).
* 93% overall alignment rate -- fantastic.

### Gene Counts (featureCounts)
And finally, we counted the genes with featureCounts. We instruct it to not count anything that  aligned to more than known gene, and to allow cases where only one read mapped. Result?
```
Status  aligned.bam
Assigned        33501959
```
**33M reads/pairs counted** for this example over the 46,983 annotated mouse genes in Ensembl. Perfect.

## Screen 2 QC and alignment stats
For screen2, the 4 sample, ribotagged, dissected experiment, the QC alerted us to the presence of unknown (to our trimming software) adapter/primer sequences and PCR-related sequences (oligo-dT) that likely originate from the particular kits used to amplify from low-input RNA. However, since we quantify only what successfully and uniquely aligns to the mouse genome, we did just that, and have left the precise origin of the artifact to be determined.

### QC (fastQC)
* http://bioinformatics.unr.edu/share/jr/gould_qc/read_2.P.trim_fastqc.html
* Basic info: ~41M reads, 50 letters long per half, Illumina 1.5 (phred64) file encoding
* Many more tests fail or warn, stemming from adapter fragments and polyTs (see _Overrepresented Sequences_ below)

#### Overrepresented Sequences
| % | sequence |
| ------| -----------|
|3.0%|**AAGCAGTGGTATCAACGCAGAGTAC**TTTTTTTTTTTTTTTTTTTTTTTTT|
|0.6%|TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT|
|0.5%|**AAGCAGTGGTATCAACGCAGAGTAC**ATGGGGAGGCATTCAGGCAGCGAGA|
|0.1%|**AAGCAGTGGTATCAACGCAGAGTAC**ATGGGAGGCATTCAGGCAGCGAGAG|
|0.1%|**AAGCAGTGGTATCAACGCAGAGTAA**ACACATCTTGCAGAGGAAGGAGGCT|
* We ran the 3rd & 5th item on that list through [vecscreen](http://www.ncbi.nlm.nih.gov/tools/vecscreen/)
* **Result:** 24-25 first nts are from some known adapter that was presumably used somewhere used in the preps.
![vecscreen_results](http://bioinformatics.unr.edu/share/jr/gould_qc/screen2-oligo.png)
why the 26-50 bases are yellow "suspects"
http://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Suspect
* However, since we align to the mouse genome, and count unique alignments, non-mouse sequences will fail to align, and the only expected harm is a reduced success rate.

### Alignment (HISAT2)
After trimming and re-running fastqc on the trimmed reads, we aligned all 41M pairs to the mouse genome (and all known splice sites) with HISAT2, and the result summary looked great:
```
25698754 reads; of these:
  25698754 (100.00%) were paired; of these:
    5744582 (22.35%) aligned concordantly 0 times
    13187114 (51.31%) aligned concordantly exactly 1 time
    6767058 (26.33%) aligned concordantly >1 times
    ----
    5744582 pairs aligned concordantly 0 times; of these:
      1161558 (20.22%) aligned discordantly 1 time
    ----
    4583024 pairs aligned 0 times concordantly or discordantly; of these:
      9166048 mates make up the pairs; of these:
        4527894 (49.40%) aligned 0 times
        2852305 (31.12%) aligned exactly 1 time
        1785849 (19.48%) aligned >1 times
91.19% overall alignment rate
```
* 13 M pairs both aligned exactly once.
* 3 M solo reads aligned once (partner read had zero alignments)
* 91% overall alignment rate.

### Gene Counts (featureCounts)
Finally, we counted the genes with featureCounts. We instruct it to not count anything that  aligned to more than known gene, and to allow cases where only one read mapped. Result?
```
Status  aligned.bam
Assigned        16391870
```
**16.4 M pairs or solo reads were counted** for this example over all 46,983 annotated mouse genes in Ensembl. _Workable._
