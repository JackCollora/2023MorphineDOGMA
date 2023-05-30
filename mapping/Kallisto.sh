#point of this is to map the reads much faster using kallisto psuedoalignments (which are as accurate or more so than traditional alignments)

#first follow the kite install instructions

#then generate index
#start with a file the same as for cell ranger 
#mine is "dogma_antibody_list_with_barcode_Yulong.csv"

#in python convert that to a name,barcode pair using this python code from the kite github
import csv 

def get_tags(input_ref, output_ref):
    with open(input_ref, mode='r') as csv_in:
        with open(output_ref, mode='w', newline='') as csv_out:
            csv_reader = csv.reader(csv_in)
            csv_writer = csv.writer(csv_out, delimiter=',')
            csv_writer.writerow(['Feature Barcode names', 'Feature Barcode sequences'])
            next(csv_reader)
            for row in csv_reader:
                csv_writer.writerow([row[1].strip(), row[4].strip()])            
    return

get_tags("dogma_antibody_list_with_barcode_Yulong.csv", "FeatureBarcodes.csv")

#convert to a fasta including all one basepair mismatches 
python ~/project/bigbin/kite/featuremap/featuremap.py FeatureBarcodes.csv --header

#make Kallisto index 

kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa

#grab the other reference file you need - barcode whitelist 

zcat /vast/palmer/apps/avx2/software/CellRanger-ARC/2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz > 10xmulti_whitelist.txt

#from here you do the mapping, filtering, bc correcting, etc. Below are snippits from the kite repository, adapted to just one sample (JC4)
kallisto bus -i FeaturesMismatch.idx -o CITE_res/JC4/ -x 10xv3 -t 10 ~/CITE_reads/"C-JC4_"*R1* ~/CITE_reads/"C-JC4_"*R3*
bustools correct -w CITEref/10xmulti_whitelist CITE_res/JC4/output.bus -o CITE_res/JC4/output_corrected.bus
bustools sort -t 10 -o CITE_res/JC4/output_sorted.bus CITE_res/JC4/output_corrected.bus
mkdir CITE_res/JC4/featurecounts
bustools count -o CITE_res/JC4/featurecounts/featurecounts --genecounts -g CITEref/FeaturesMismatch.t2g -e CITE_res/JC4/matrix.ec -t CITE_res/JC4/transcripts.txt CITE_res/JC4/output_sorted.bus


#final job file for all of them
#JC 7 and 11 have an _ instead of a - for some reason, manually corrected afterwards in the job file
for i in JC4 JC5 JC6 JC7 JC8 JC9 JC10 JC11 JC12; do echo module load miniconda";" conda activate kite";" mkdir CITE_res/$i ";" kallisto bus -i CITEref/FeaturesMismatch.idx -o CITE_res/$i/ -x 10xv3 -t 10 CITE_reads/"C-"$i"_"*R1* CITE_reads/"C-"$i"_"*R3* ";" bustools correct -w CITEref/10xmulti_whitelist.txt CITE_res/$i/output.bus -o CITE_res/$i/output_corrected.bus";" bustools sort -t 10 -o CITE_res/$i/output_sorted.bus CITE_res/$i/output_corrected.bus";" mkdir CITE_res/$i/featurecounts ";" bustools count -o CITE_res/$i/featurecounts/featurecounts --genecounts -g CITEref/FeaturesMismatch.t2g -e CITE_res/$i/matrix.ec -t CITE_res/$i/transcripts.txt CITE_res/$i/output_sorted.bus; done


#after this runs, then its just a matrix file like any other 10X. 
for i in $(seq 4 12); do cp -r ~/scratch/10X/CITE_res/JC$i/featurecounts JC$i/CITEmatrix; done
