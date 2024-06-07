#!/bin/bash

##################################################################
##################################################################
#
# 	19 May 2024 revision:
#	-- remove step in combining assembled and unassembled contig
#	-- analyze only assembled contig. but keep copy of unassembled contig	
#
#	Modidifications from version include:
#	-- Add Step 5: Classify contigs by local blastn
#	-- Do not delete / keep prinseq-trimmed fastq
#	-- all tools and databases should be specified here first by user; dependent on machine used; here is for Alpha2
#
#	Name of program: canu_metagenomics_v1.sh
# 	Created by: Majo Galarion
# 	Date created: 20 February 2024
# 	Steps:
#	Filters reads with <200 bp and <9 mean Q-score
#	Corrects and assembles read with Canu
#	Classifies contigs against protein db -- diamond blastx
#	Classifies contigs with kraken2
#	
#	
#	This script can be used on demultiplexed ONT metagenomics datasets
#	
#
##################################################################
##################################################################

# Always run this script inside directory of your input files!

### Important dependencies:
# Define tools and database and how you call them in your current system
# Differs with every machine / system
#
prinseq="prinseq-lite.pl"
canu="canu"
diamond="diamond"
diamdb="/db/diamond/ViralRefSeqProtein.dmnd2"
krona="/home3/2509094g/KronaTools-2.8.1/bin/bin"
kraken="/software/kraken2-v2.1.1/kraken2"
krakendb="/home3/2509094g/kraken-DB/minikraken2_v2/minikraken2_v2_8GB_201904_UPDATE"
blastn="blastn"
blastdb="/db/blast_v5/nt"


usage() {
	echo "";
	echo "ONT metagenomics analysis";
	echo "";
	echo "Usage $0 -i reads.fastq" 1>&2;
	echo "";
	echo "Optional parameters:";
	echo "-t: number of threads (default: 8)";
	echo "-m: minimum read length (default: 200)";
	echo "-x: maximum read length (default: none)";
	echo "-d: minimum read depth (default: 20)";
	echo "-q: minimum read Q-score (default: 9)";
	echo "-s: estimated genome size (default: 10k)";
	echo "";
	exit 1;
}

if [ $# -eq 0 ]; then usage; exit 1; fi


### These are the default paramters
threads="8";
minLength=200;
maxLength="none";
minDepth=20;
minQual=9;
genomeSize="10k";
blast_max_target="1";



### Read these arguments when specified in command
while getopts ":i:t:m:x:d:q:s:" opt
do
	case "${opt}" in
		i) fq=${OPTARG} ;;
		t) threads=${OPTARG} ;;
		m) minLength=${OPTARG} ;;
		x) maxLength=${OPTARG} ;;
		d) minDepth=${OPTARG} ;;
		q) minQual=${OPTARG} ;;
		s) genomeSize=${OPTARG} ;;
		b) blast_max_target=${OPTARG} ;;
		?) echo "Option -${OPTARG} requires an argument" >&2
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

### Define name of sample to be used as prefix
name=${fq%.fastq.gz}
name=${name%.fq.gz}
name=${name%.fq}
name=${name%.fastq}

### Create output directory
mkdir meta_${name}

### Output arguments specified to terminal stdout
echo "";
echo -e "\033[31m ONT metagenomics analysis\033[0m";
echo -e "\033[31m by Majo Galarion\033[0m";
echo "";
echo -e "Input file: \033[33m${fq}\033[0m";
echo -e "Number of threads: \033[33m${threads}\033[0m";
echo -e "Min read length: \033[33m${minLength}\033[0m";
echo -e "Max read length: \033[33m${maxLength}\033[0m";
echo -e "Minimum depth: \033[33m${minDepth}\033[0m"
echo -e "Minimum read Q-score: \033[33m${minQual}\033[0m"
echo -e "Estimated genome size: \033[33m${genomeSize}\033[0m"
echo -e "Blast maximum target: \033[33m${blast_max_target}\033[0m"

##################################################################



### Step 1: Size selection with prinseq-lite
echo "";
echo -e "\033[32m Step 1/5: Size selection with prinseq-lite...\033[0m"
echo "";

	${prinseq} -fastq ${fq} -min_len ${minLength} -min_qual_mean ${minQual} -out_format "3" -out_good meta_${name}/${name}.prinseq
	rm -f *_prinseq_bad_*

	# get read length distribution on raw reads before and after porechop and prinseq
	cat ${fq} | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > meta_${name}/${fq}_read-len-dist.csv
	cat meta_${name}/${name}.prinseq.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > meta_${name}/${name}.prinseq_read-len-dist.csv


### Step 2: Run canu to correct reads and generate assembly
echo "";
echo -e "\033[32m Step 2/5: Correct reads and generate assembly using Canu...\033[0m"
echo "";

	# define overlap length for assembly; should be at least half of min length
	overlap="${minLength}/4"
	${canu} -nanopore-corrected meta_${name}/${name}.prinseq.fastq -d meta_${name}/canu_correct_assemble -p ${name}.canu useGrid=false minReadLength=${minLength} minOverlapLength=${overlap} genomesize=${genomeSize} corMinCoverage=0

	# combine unassemled and assembled contigs
	# keep original contigs in canu output dir
	#cat meta_${name}/canu_correct_assemble/${name}.canu.contigs.fasta meta_${name}/canu_correct_assemble/${name}.canu.unassembled.fasta > meta_${name}/${name}.canu.all-contigs.fasta
	#
	# copy assembled contig for analysis
	cp meta_${name}/canu_correct_assemble/${name}.canu.contigs.fasta meta_${name}/${name}.canu.all-contigs.fasta	

# create a summary table for all contigs generated by canu
grep ">" meta_${name}/${name}.canu.all-contigs.fasta | while read -r line;
do
	contig=$(echo ${line} | awk '{print $1}' | cut -c2-)
	len=$(echo ${line} | awk '{print $2}' | cut -c5-)
	reads=$(echo ${line} | awk '{print $3}' | cut -c7-)
	echo ${name} ${contig} ${len} ${reads} | tr " " "\t" > meta_${name}/${contig}.stats.temp
done

cat meta_${name}/*.stats.temp > meta_${name}/${name}.canu.contiglist.txt
sed -i '1i Sample\tContig_ID\tLength\tReads' meta_${name}/${name}.canu.contiglist.txt
rm meta_${name}/*.stats.temp


# simplify contig headers
cat meta_${name}/${name}.canu.all-contigs.fasta | awk '/^>/{print $1}!/>/' > meta_${name}/${name}.canu.contigs.fasta
rm meta_${name}/${name}.canu.all-contigs.fasta


### Step 3: Classify contigs using diamond blastx
echo "";
echo -e "\033[32m Step 3/5: Classify contigs with diamond blastx...\033[0m"
echo "";

	mkdir meta_${name}/diam_temp
	${diamond} blastx -d ${diamdb} -f 6 -o meta_${name}/${name}.diamond.viral.txt -p ${threads} -q meta_${name}/${name}.canu.contigs.fasta -t meta_${name}/diam_temp

	# import taxonomy to krona
	${krona}/ktImportBLAST meta_${name}/${name}.diamond.viral.txt -p -k -o meta_${name}/${name}.krona.viral.html


### Step 4: Classify contigs using Kraken2
echo "";
echo -e "\033[32m Step 4/5: Classify contigs with kraken2...\033[0m"
echo "";

${kraken} --db ${krakendb} --threads ${threads} --use-names --output meta_${name}/${name}.kraken.output --report meta_${name}/${name}.kraken.report meta_${name}/${name}.canu.contigs.fasta

# get kraken output summary for plotting
grep "tig0" meta_${name}/${name}.kraken.output | while read -r line;
do
	contig=$(echo ${line} | awk '{print $2}')
	taxid=$(echo ${line} | awk '{print $3,$4,$5,$6,$7}' | sed 's/(.*//') 
	echo ${name} ${contig} ${taxid} | tr " " "\t" > meta_${name}/${contig}.kraken.taxid.temp
done

cat meta_${name}/*.kraken.taxid.temp > meta_${name}/${name}.kraken.taxid.txt
sed -i '1i Sample\tContig_ID\tTaxa_classification' meta_${name}/${name}.kraken.taxid.txt
rm meta_${name}/*.kraken.taxid.temp


### Step 5: Classify contigs using blastn
echo "";
echo -e "\033[32m Step 5/5: Classify contigs with blastn...\033[0m"
echo "";

${blastn} -db ${blastdb} -query meta_${name}/${name}.canu.contigs.fasta -out meta_${name}/${name}.blastn.report.txt -max_target_seqs ${blast_max_target} -num_threads ${threads} -outfmt '6 qseqid sacc qstart qend sstart send evalue bitscore length pident nident mismatch sscinames scomnames sblastnames sskingdoms'

sed -i '1i Query_seq_id\tSubject_accession\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tBistscore\tLength\tPerc_identical_match\tNum_identical_match\tMismatch	Subject_SciName\tCommon_name\tBlast_name\tSuper_kingdom' meta_${name}/${name}.blastn.report.txt


# Don't remove prinseq-trimmed fastq as can be used for ref-based assembly
#rm meta_${name}/${name}.prinseq.fastq

rm -r meta_${name}/diam_temp
rm  meta_${name}/${name}.kraken.output

echo "";
echo -e "\033[31m Metagenomics nalysis done!!!\033[0m";
echo "";