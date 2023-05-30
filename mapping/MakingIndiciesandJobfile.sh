#running cellranger

#generating libraries.csv 
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo "fastqs,sample,library_type" >$i"_library.csv"; echo ~/palmer_scratch/10X/RNA_reads/,$i"_MAH",Gene Expression >>$i"_library.csv"; echo ~/palmer_scratch/10X/ATAC_reads/,"A-"$i"_MAH",Chromatin Accessibility >>$i"_library.csv"; done 

#making a job file 

for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load CellRanger-ARC";" cellranger-arc count --id=$i --reference=/home/jac369/project/refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --libraries=librarie/$i"_library.csv" --localcores=20 --localmem=120; done

#make it into a dSQ script
dSQ --jobfile cellrangerjobs.txt -c 20 --mem 125g -p day -t 1- -N 1 

#submit
sbatch dsq-cellrangerjobs-2023-04-26.sh

#making indicies for ATAC and RNA

module load Bowtie2; bowtie2-build --threads 10 HXB2.fasta bowtie_index/HXB2
module load Bowtie2; bowtie2-build --threads 10 cladeB.fasta bowtie_index/cladeB
module load Bowtie2; bowtie2-build --threads 10 morphine_autologous.fasta bowtie_index/morphine_autologous
module load STAR; STAR --runThreadN 9 --runMode genomeGenerate --genomeDir STAR_index/HXB2 --genomeFastaFiles HXB2.fasta --genomeSAindexNbases 5 --outTmpDir crap2
module load STAR; STAR --runThreadN 9 --runMode genomeGenerate --genomeDir STAR_index/CladeB --genomeFastaFiles cladeB.fasta --genomeSAindexNbases 12 --outTmpDir crap3
module load STAR; STAR --runThreadN 9 --runMode genomeGenerate --genomeDir STAR_index/morphine_autologous --genomeFastaFiles morphine_autologous.fasta --genomeSAindexNbases 7 --outTmpDir crap1

####second iteration
#then we actually want to run the STAR command, just to get it going 

for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load STAR";" STAR --runThreadN 10 --runMode alignReads --genomeDir ~/project/refs/STAR_index/HXB2 --readFilesIn RNA_reads/$i"_R3_Cat.fastq.gz" --outFilterMultimapNmax 1000 --outFileNamePrefix HIV_RNA/$i"_HXB2" --readFilesCommand zcat; done >RNA_mapping.txt
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load STAR";" STAR --runThreadN 10 --runMode alignReads --genomeDir ~/project/refs/STAR_index/CladeB --readFilesIn RNA_reads/$i"_R3_Cat.fastq.gz" --outFilterMultimapNmax 10000 --outFileNamePrefix HIV_RNA/$i"_cladeB" --readFilesCommand zcat; done >>RNA_mapping.txt
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load STAR";" STAR --runThreadN 10 --runMode alignReads --genomeDir ~/project/refs/STAR_index/morphine_autologous --readFilesIn RNA_reads/$i"_R3_Cat.fastq.gz" --outFilterMultimapNmax 1000 --outFileNamePrefix HIV_RNA/$i"_autologous" --readFilesCommand zcat; done >>RNA_mapping.txt

#mapping the ATAC with bowtie2

for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load Bowtie2 SAMtools";" bowtie2 -p 10 --very-sensitive --no-unal -a -x ~/project/refs/bowtie_index/cladeB -1  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R1* -2  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R3* - "|" samtools view -u - "|" samtools sort - -o HIV_ATAC/$i"CladeB.bam"; done >ATACjobs.txt
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load Bowtie2 SAMtools";" bowtie2 -p 10 --very-sensitive --no-unal -a -x ~/project/refs/bowtie_index/HXB2 -1  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R1* -2  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R3* - "|" samtools view -u - "| "samtools sort - -o HIV_ATAC/$i"HXB2.bam"; done >>ATACjobs.txt
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load Bowtie2 SAMtools";" bowtie2 -p 10 --very-sensitive --no-unal -a -x ~/project/refs/bowtie_index/morphine_autologous -1  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R1* -2  ~/scratch/10X/ATAC_reads/"A-"$i"_MAH"*R3* - "|" samtools view -u - "|" samtools sort - -o HIV_ATAC/$i"autologous.bam"; done >>ATACjobs.txt


