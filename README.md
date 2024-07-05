# DENV metagenomics
This script is a custom metagenomics pipeline that performs de novo assembly and classification on uncompressed raw FASTQ files<br/>
by **M Galarion**
<br/>
<br/>
## IMPORTANT DEPENDENCIES
Make sure the following tools are installed:<br/>
&emsp;**PRINSEQ-lite**  https://sourceforge.net/projects/prinseq/files/standalone/<br/>
&emsp;**canu**  https://github.com/marbl/canu<br/>
&emsp;**kraken2**  https://github.com/DerrickWood/kraken2<br/>
&emsp;**kraken database**<br/>
&emsp;**local BLASTn** https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/<br/>
&emsp;**BLASTn database**<br/>
<br/>
Once these tools are installed, you need to explicitly define the command usage for each tool and specify the location of the databases <br/>
To do this, open the latest version of the bash script in a text editor <br/>
Starting from around line 40, edit the line that corresponds to the tool definition (inside the quotation marks) <br/>
 <br/>
For example:
```

### Define tools and how you call them in your current system

prinseq = "prinsesq-lite.pl"
canu = "canu"
kraken = "/software/kraken2-v2.1.1/kraken2"
krakendb = "/home/user/kraken-DB/minikraken2_v2/minikraken2_v2_8GB_201904_UPDATE"    # directory of kraken2 database
blastn = "blastn"
blastdb = "/db/blast_v5/nt"    # directory of blast nt database

```
The tool definition inside the quotation marks should correspond to how you call them in your current system<br/>
Save the bash file and close<br/>

<br/>

### The script is an executable file written in bash and performs the following steps:<br/>
1. Trims the raw FASTQ files based on user-specified minimum length and minimum read Q-score <br/>
2. Read correction and de novo assembly using canu <br/>
3. K-mer based classification of assembled contigs using Kraken2 and specified kraken2 database <br/>
4. Pairwise alignment of assembled contigs using BLASTn and specified blast database <br/>
5. Report generation <br/>
<br/>

### Command line usage:
```

bash canu_metagenomics_v2.1.sh -i raw_reads.fastq

```
**NOTE:** Always run the script inside the directory where your -i files are located.

<br/>

### To print options and default values:
```

bash canu_metagenomics_v2.1.sh

```
```
Optional parameters:
-t: number of threads (default: 8)
-m: minimum read length (default: 200)
-x: maximum read length (default: none)
-d: minimum read depth (default: 20)
-q: minimum read Q-score (default: 9)
-s: estimated genome size (default: 10k)
```
