#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

sample=$1 #name of sample directory
Read1=$(findFiles.sh "${sample}" | awk 'NR==1{print $1}') #name of read1
Read2=$(findFiles.sh "${sample}" | awk 'NR==2{print $1}') #name of read2

################################################################
# Utilities
################################################################

# Record hostname
function getHostname {
	cd "${Data}"/"${sample}"
	hostname > host.txt
}

# Function timer
function funTime {
	local func=$1
	echo "" >> "${Data}"/"${sample}"/Runtime."${sample}".log
	echo "${func} started" >> "${Data}"/"${sample}"/Runtime."${sample}".log
	{ time "${func}"; } >> "${Data}"/"${sample}"/Runtime."${sample}".log 2>&1
	echo "${func} done" >> "${Data}"/"${sample}"/Runtime."${sample}".log
}

################################################################
# Quality Control of experimental data
################################################################

# FastQC
function runfastq {
	cd "${Data}"/"${sample}"/
	mkdir fastqc"${runNumber}"."${sample}"
	fastqc -t $PROCS -o fastqc"${runNumber}"."${sample}" *."$fastq".gz
	cd "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/
	unzip \*.zip
}

#Scraping Illumina Version Number to get phred number
function getPhred {
	cd "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/*/
	IllumVN=$(grep -E "Encoding" fastqc_data.txt | sed 's/[^0-9.]*//g')
	echo $(grep -E "Encoding" fastqc_data.txt)
	if [ `bc <<< "$IllumVN >= 1.3 && $IllumVN < 1.8"` -eq 1 ]; then
		phred="64"
	else
		phred="33"
	fi
}

#Scraping Overrepresented Sequences
function findOverrepSeq {
	cd "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/
	for directory in $(ls -d */); do
		cd "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/"$directory"/
		awk '/Overrepresented sequences/,/>>END_MODULE/' fastqc_data.txt \
			| tail -n +2 | head -n -1 > Overrepresented_sequences.txt
		if [ $(wc -l < Overrepresented_sequences.txt) -eq 0 ]; then
			echo "There are no over represented sequences in "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/"$directory""
		else
			echo "There are $overseq over represented sequences in "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/"$directory"\
				; see Overrepresented_sequences.txt in "${Data}"/"${sample}"/fastqc"${runNumber}"."${sample}"/"$directory""
			#TODO exit?
		fi
	done
}

# Trimmomatic
function runTrim {
	cd "${Data}"/"${sample}"/
	java -jar $RNASEQDIR/Trimmomatic/trimmomatic-0.35.jar \
		PE -threads $PROCS -phred"${phred}" \
		"${Read1}" "${Read2}" \
		read1.P.trim."$fastq".gz read1.U.trim."$fastq".gz \
		read2.P.trim."$fastq".gz read2.U.trim."$fastq".gz \
		ILLUMINACLIP:$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25
}

# Running fastqc, getting overrepresented sequences, TODO doing something, running trimmomatic
function qcData {
	runNumber=1
	funTime runfastq $runNumber
	if [ $runNumber -eq 1 ]; then
		funTime getPhred $runNumber
	else
		:
	fi
	funTime findOverrepSeq $runNumber
	#TODO Do something to help trimmomatic
	funTime runTrim
}

#TODO wrapper to automate qcData

################################################################
# Pipeline
################################################################

# SeqTK sub sample
function runSeqtk {
	seqtk sample -s100 read1.P.trim."$fastq".gz 10000 | seqtk seq -A - > sampled.read1.fa
	seqtk sample -s100 read2.P.trim."$fastq".gz 10000 | seqtk seq -A - > sampled.read2.fa
}

# BLAST the sub samples
function runBlastn {
	blastn -query sampled.read1.fa \
		-db ""${Ref}"/"$basename".cdna.all" \
		-out sampled.read1_vscdna.out \
		-task blastn-short -outfmt '6 std sstrand' \
		-max_target_seqs 1 -num_threads $PROCS
	blastn -query sampled.read2.fa \
		-db ""${Ref}"/"$basename".cdna.all" \
		-out sampled.read2_vscdna.out \
		-task blastn-short -outfmt '6 std sstrand' \
		-max_target_seqs 1 -num_threads $PROCS
}

# Stranded classifier.  Should only continue if True.  NOT EVALUATED.
function findStranded {
	stranded_classifier.py -1 sampled.read1_vscdna.out -2 sampled.read2_vscdna.out
	#TODO do something if stranded or not
}

# Hisat2
function runHisat {
	# stranded classifier stuff that gets done
	hisat2 -k 5 -p $PROCS --dta --phred"$phred" \
	--known-splicesite-infile "${Ref}"/splice_sites.txt \
	-x "${Ref}"/"$basename" \
	-1 read1.P.trim."$fastq".gz -2 read2.P.trim."$fastq".gz \
	-S aligned."${sample}".sam
}

# SAM -> BAM
function runCompression {
	samtools view -bT "${Ref}"/"${Genome}" \
	-@$PROCS aligned."${sample}".sam -o aligned."${sample}".bam
	rm aligned."${sample}".sam
}

# To counts
function runFeatureCounts {
	featureCounts -T $PROCS -p -C --primary --ignoreDup -t exon -g gene_id \
	-a "${Ref}"/"${Gtf}" -o aligned."${sample}".counts aligned."${sample}".bam
}

# Puts geneid and aligned.bam columns into a separate file
function getKarensColumns {
	tail -n +2 aligned."${sample}".counts | awk '{printf ("%5s\t%s\t%s\n", $1, $6, $7)}' > aligned."${sample}".counts.three
}

function getAlignedColumn {
	tail -n +2 aligned."${sample}".counts | awk '{print $7}' > aligned."${sample}".counts.one
}


function runPipeline {
	cd "${Data}"/"${sample}"/
	funTime runSeqtk
	funTime runBlastn
	funTime findStranded
	funTime runHisat
	funTime runCompression
	funTime runFeatureCounts
	funTime getKarensColumns
	funTime getAlignedColumn
}

function main {
	echo "Runtime Log" >> "${Data}"/"${sample}"/Runtime."${sample}".log
	echo "----------------------------------------" >> "${Data}"/"${sample}"/Runtime."${sample}".log
	funTime getHostname
	#TODO make wrapper for qcData and call that instead
	#	wrapper will automate iterations of qcData
	funTime qcData
	funTime runPipeline
}

> "${Data}"/"${sample}"/Runtime."${sample}".log
funTime main
