#!/bin/bash

Usage="QCofRef.sh [-h]

where:
    -h          Show this screen"

while getopts ':h' option; do                           
    case "${option}" in                                 
        h) echo "${Usage}"; exit;;                      
        ?)                                              
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;                     
    esac                                                
done                                                    
shift $(( OPTIND -1 ));                                 

if [ -f "${Input_Field}" ]; then
    source "${Input_Field}"
else
    echo "NEED INPUT FILE"
    exit
fi

################################################################
# Quality Control of reference data
################################################################

function qcReport {
    cd "${Ref}"
    
    #Header
    echo "Reference Report" > Reference_Report.txt
    echo "----------------------------------------------" >> Reference_Report.txt
    
    #Getting Chromosome names and count in Genome and GTF
    grep -E '^>' "${Genome}" | awk '{print $1}' | cut -sf 2 -d '>' | sort > gen_chrom_names.txt
    gen_chrom_count=$(wc -l gen_chrom_names.txt | cut -sf 1 -d " ")
    
    cut -sf 1 *.gtf | uniq | sort > gtf_chrom_names.txt
    gtf_chrom_count=$(wc -l gtf_chrom_names.txt | cut -sf 1 -d " ")
    
    #Chromosome names sample
    echo "Chromosome names:" >> Reference_Report.txt
    echo "GENOME    GTF" >> Reference_Report.txt
    paste gen_chrom_names.txt gtf_chrom_names.txt >> Reference_Report.txt
    
    echo "" >> Reference_Report.txt
    echo "Total number of Chromosomes:" >> Reference_Report.txt
    echo "$gen_chrom_count  $gtf_chrom_count" >> Reference_Report.txt
    
    #diff of Chromosome names
    diff gen_chrom_names.txt gtf_chrom_names.txt > chrom_diff.txt
    chromDiff=$(wc -l chrom_diff.txt | cut -sf 1 -d " ")
    if [ "$chromDiff" -eq 0 ]; then
        :
    else
        echo "CHROMOSOME NAMES NEED HELP; see chrom_diff.txt" >> Reference_Report.txt
    fi
    
    #Getting genes in cDNA and GTF
    grep -Eo '^>\S+' *cdna* | cut -sf 2 -d '>' | cut -sf 1 -d '.' \
        | uniq | sort | uniq > genes_in_cDNA.txt
    genesInCdna=$(wc -l genes_in_cDNA.txt | cut -sf 1 -d " ")
    cut -sf 9 *.gtf | cut -sf 1 -d ';' | cut -sf 2 -d '"' \
        | uniq | sort | uniq > genes_in_gtf.txt
    genesInGTF=$(wc -l genes_in_gtf.txt | cut -sf 1 -d " ")
    
    
    #Gene names sample
    echo "" >> Reference_Report.txt
    echo "Sample of gene names:" >> Reference_Report.txt
    echo "cDNA      GTF" >> Reference_Report.txt
    (paste genes_in_cDNA.txt genes_in_gtf.txt \
        | head -n 5 >> Reference_Report.txt) 2> /dev/null
    
    #Getting number of transcript names in cDNA
    echo "" >> Reference_Report.txt
    echo "Total number of gene transcripts in cDNA:" >> Reference_Report.txt
    grep -Ec '^>' *cdna* >> Reference_Report.txt
    
    
    echo "" >> Reference_Report.txt
    echo "Total number of unique gene names:" >> Reference_Report.txt
    echo "$genesInCdna      $genesInGTF" >> Reference_Report.txt
    
    #Finding gene names in cDNA but not in GTF
    diff genes_in_cDNA.txt genes_in_gtf.txt > gene_diff.txt
    geneDiff=$(diff genes_in_cDNA.txt genes_in_gtf.txt | grep '<' | wc -l)
    
    echo "" >> Reference_Report.txt
    echo "Number of gene names in cDNA not in GTF: $geneDiff" >> Reference_Report.txt
}


cd "${Ref}"
qcReport
cat Reference_Report.txt
