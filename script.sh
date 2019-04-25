#1- Data Download
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra
fastq-dump SRR8797509.sra

#2- Prepare the data
#a- downsample the fastq file from 23 million reads to only 5 million reads
seqkit head -n 5000000 SRR8797509.fastq > sample_5_mil.fastq
#b- split the fastq file to 5 separate files
seqkit split sample_5_mil.fastq -p 5 -O unshuffled
#c- shuffle the reads and split them again into 5 files
seqkit shuffle sample_5_mil.fastq | seqkit split -p 5 -O shuffled

#3- quality control 
fastqc -f fastq -noextract shuffled/stdin.part_001.fastq unshuffled/sample_5_mil.part_001.fastq -o qc


#4-trimming 
mkdir -p trimmed trimmed/mild trimmed/aggressive
#mild trimming
for f1 in unshuffled/*.fastq
do
 f2=$(echo "$f1"| sed 's,unshuffled,trimmed/mild,g')
 trimmomatic SE -phred33 $f1 ${f2%fastq}"trimmed.fq.gz" ILLUMINACLIP:$adap/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done
#aggresive trimming
for f1 in shuffled/*.fastq
do
 f2=$(echo "$f1"| sed 's,shuffled,trimmed/aggressive,g')
 trimmomatic SE -phred33 $f1 ${f2%fastq}"trimmed.fq.gz" ILLUMINACLIP:$adap/TruSeq2-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50
done

#5-alignment
#bwa
mkdir bwa_align && cd bwa_align

bwa index -a bwtsw GRCh38.p12.genome.fa

for input in ../unshuffled/*.fastq
do
 output=$(echo "$input" | sed 's,fastq,sam,g')
 bwa mem GRCh38.p12.genome.fa $input > $output
done
#hisat
mkdir hisat_align && cd hisat_align

hisat2_extract_splice_sites.py ../gencode.v30.annotation.gtf > splicesites.tsv
hisat2_extract_exons.py ../gencode.v30.annotation.gtf > exons.tsv
hisat2-build --ss splicesites.tsv --exon exons.tsv GRCh38.p12.genome.fa GRCh38.p12.genome

