#!/bin/bash
#
## This script is intended for making consensus sequences
readarray -t barcodes < barcodes.txt
#for folder in barcode*
for folder in "${barcodes[@]}"
do
	cd $folder
	cp ../consensus_seq_and_score.denovo2.R .
	Rscript consensus_seq_and_score.denovo2.R Match1.fasta minimap_cluster1.maf.cons.data
	rm consensus_seq_and_score.denovo2.R
	cd ..
done
mkdir results2
cp ./barcode*/*.barcode*.*.cons.denovo.python.fasta results2
cd results2
cat *.fasta > all_seqs.fasta
### blast
blastn -query all_seqs.fasta -db ~/blastdb/nt -max_target_seqs 1 -num_threads 16 -outfmt "6 qseqid sacc pident length mismatch evalue bitscore staxids sscinames stitle" -out blast.out
