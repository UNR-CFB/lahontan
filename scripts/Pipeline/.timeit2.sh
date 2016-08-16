#!/bin/bash

if [ -f /home/alberton/Pipeline/INPUT ]; then
	source /home/alberton/Pipeline/INPUT
else
	echo "NEED INPUT FILE"
	exit
fi

Read1=$1 #name of read1
Read2=$2 #name of read2

# Record hostname
function getHostname {
	cd "${Data}"/
	hostname > host.txt
}

################################################################
# Utilities
################################################################

function funTime {
	local func=$1
	echo "${func}" >> "${Data}"/Runtime.log
	{ time "${func}"; } >> "${Data}"/Runtime.log 2>&1
}

################################################################
# Quality Control of experimental data
################################################################

# FastQC
function runfastq {
	#local runNumber=$1
	cd "${Data}"/
	mkdir fastqc"${runNumber}"
	fastqc -t $PROCS -o fastqc"${runNumber}" *."$fastq".gz
	cd "${Data}"/fastqc"${runNumber}"/
	unzip \*.zip
}

#Scraping Illumina Version Number to get phred number
function getPhred {
	#local runNumber=$1
	cd "${Data}"/fastqc"${runNumber}"/*/
	IllumVN=$(grep -E "Encoding" fastqc_data.txt | sed 's/[^0-9.]*//g')
	echo $(grep -E "Encoding" fastqc_data.txt)
	if [ `bc <<< "$IllumVN >= 1.3 && $IllumVN < 1.8"` -eq 1 ]; then
		phred="64"
		echo 64 > phrednum
	else
		phred="33"
		echo 33 > phrednum
	fi
	echo "################### this is phred: ${phred}"
}

#Scraping Overrepresented Sequences
function findOverrepSeq {
	#local runNumber=$1
	cd "${Data}"/fastqc"${runNumber}"/
	for directory in $(ls -d */); do
		cd "${Data}"/fastqc"${runNumber}"/"$directory"/
		awk '/Overrepresented sequences/,/>>END_MODULE/' fastqc_data.txt \
			| tail -n +2 | head -n -1 > Overrepresented_sequences.txt
		if [ $(wc -l < Overrepresented_sequences.txt) -eq 0 ]; then
			echo "There are no over represented sequences in "${Data}"/fastqc"${runNumber}"/"$directory"/"
		else
			echo "There are $overseq over represented sequences in "${Data}"/fastqc"${runNumber}"/"$directory"/\
				; see Overrepresented_sequences.txt in "${Data}"/fastqc"${runNumber}"/"$directory"/"
			#TODO exit?
		fi
	done
	echo "################### this is phred: ${phred}"
}

# Trimmomatic
function runTrim {
	cd "${Data}"/
	echo "################### this is phred: ${phred}"
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
	#local runNumber=$1
	runNumber=1
	funTime runfastq $runNumber
	if [ $runNumber -eq 1 ]; then
		funTime getPhred $runNumber
	else
		:
	fi
	funTime findOverrepSeq $runNumber
	echo "################### this is phred: ${phred}"
	#TODO Do something to help trimmomatic
	funTime runTrim
}

#TODO wrapper to automate qcData

################################################################
# Pipeline
################################################################

# SeqTK sub sample
function runSeqtk {
	seqtk sample -s100 read1.P.trim."${fastq}".gz 10000 | seqtk seq -A - > sampled.read1.fa
	seqtk sample -s100 read2.P.trim."${fastq}".gz 10000 | seqtk seq -A - > sampled.read2.fa
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
	echo "################### this is phred: ${phred}"
	hisat2 -k 5 -p $PROCS --dta --phred"${phred}" \
		--known-splicesite-infile "${Ref}"/splice_sites.txt \
		-x "${Ref}"/"$basename" \
		-1 read1.P.trim."$fastq".gz -2 read2.P.trim."$fastq".gz \
		-S aligned.sam
}

# SAM -> BAM
function runCompression {
	samtools view -bT "${Ref}"/"${Genome}" \
	-@$PROCS aligned.sam -o aligned.bam
	rm aligned.sam
}

# To counts
function runFeatureCounts {
	featureCounts -T $PROCS -p -C --primary --ignoreDup -t exon -g gene_id \
	-a "${Ref}"/"${Gtf}" -o aligned.counts aligned.bam
}

# Puts geneid and aligned.bam columns into a separate file
function getKarensColumns {
	tail -n +2 aligned.counts | awk '{printf ("%5s\t%s\n", $1, $7)}' > aligned.counts.karen
}

function runPipeline {
	cd "${Data}"/
	funTime runSeqtk
	funTime runBlastn
	funTime findStranded
	funTime runHisat
	funTime runCompression
	funTime runFeatureCounts
	funTime getKarensColumns
}

function main {
	echo "Runtime Log" >> "${Data}"/Runtime.log
	echo "----------------------------------------" >> "${Data}"/Runtime.log
	funTime getHostname
	#TODO make wrapper for qcData and call that instead
	#	wrapper will automate iterations of qcData
	funTime qcData
	funTime runPipeline
}

echo "#" > "${Data}"/Runtime.log
funTime main
