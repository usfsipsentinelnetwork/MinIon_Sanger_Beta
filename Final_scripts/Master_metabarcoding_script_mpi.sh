#!/bin/bash
#
## This script is intended for the processing of amplicons from
## MinIon runs
#
# packages needed:
# EMBL_EBI
# MAFFT
# emboss embassy package phylip implementation (ednadist)
# python
# Databases
# mothur
# NanoFilt
# cutadapt
# R
#
# changes made 12/20/22
## changes made 12/8/22
#export MAFFT_N_THREADS_PER_PROCESS="8"
export MAFFT_N_THREADS_PER_PROCESS="64"
export MAFFT_MPIRUN="/lib64/openmpi/bin/mpirun"
# you will also need the following files in the directory
# input
# metabarcoding_consensus_files.R
# before running script conda activate ont
declare -a ARRAY_FOLDERS
ARRAY_FOLDERS[0]="18S"
ARRAY_FOLDERS[1]="BT"
ARRAY_FOLDERS[2]="COI3P"
ARRAY_FOLDERS[3]="ITS1F4"
ARRAY_FOLDERS[4]="ITS54"
ARRAY_FOLDERS[5]="ITS14"
ARRAY_FOLDERS[6]="LepCOI"
ARRAY_FOLDERS[7]="LSU"
ARRAY_FOLDERS[8]="TEF"
ARRAY_FOLDERS[9]="16S"
declare -a PRIMERS_F
declare -a PRIMERS_R
declare -a PRIMERS_FRC
declare -a PRIMERS_RRC
declare -a PRODUCT_LENGTH_MIN
declare -a PRODUCT_LENGHT_MAX
# 18S
# 965 GGCGATCAGATACCGCCCTAGTT
# 965 (revcomp) AACTAGGGCGGTATCTGATCGCC
# 1573 TACAAAGGGCAGGGACGTAAT
# 1573 (revcomp) ATTACGTCCCTGCCCTTTGTA
PRIMERS_F[0]="GGCGATCAGATACCGCCCTAGTT"
PRIMERS_FRC[0]="AACTAGGGCGGTATCTGATCGCC"
PRIMERS_R[0]="TACAAAGGGCAGGGACGTAAT"
PRIMERS_RRC[0]="ATTACGTCCCTGCCCTTTGTA"
PRODUCT_LENGTH_MIN[0]="300"
PRODUCT_LENGTH_MAX[0]="800"
# BT
# BT2A ggtaaccaaatcggtgctgctttc
# BT2A (revcomp) gaaagcagcaccgatttggttacc
# BT2B accctcagtgtagtgacccttggc
# BT2B (revcomp) gccaagggtcactacactgagggt
PRIMERS_F[1]="ggtaaccaaatcggtgctgctttc"
PRIMERS_FRC[1]="gaaagcagcaccgatttggttacc"
PRIMERS_R[1]="accctcagtgtagtgacccttggc"
PRIMERS_RRC[1]="gccaagggtcactacactgagggt"
PRODUCT_LENGTH_MIN[1]="400"
PRODUCT_LENGTH_MAX[1]="1000"
# COI3P
# A1718 GGAGGATTTGGAAATTGATTAGTTCC
# A1718 (revcomp) GGAACTAATCAATTTCCAAATCCTCC
# A2411 GCTAATCATCTAAAAACTTTAATTCCWGTWG
# A2411 (revcomp) CWACWGGAATTAAAGTTTTTAGATGATTAGC
PRIMERS_F[2]="GGAGGATTTGGAAATTGATTAGTTCC"
PRIMERS_FRC[2]="GGAACTAATCAATTTCCAAATCCTCC"
PRIMERS_R[2]="GCTAATCATCTAAAAACTTTAATTCCWGTWG"
PRIMERS_RRC[2]="CWACWGGAATTAAAGTTTTTAGATGATTAGC"
PRODUCT_LENGTH_MIN[2]="400"
PRODUCT_LENGTH_MAX[2]="1200"	
# ITS1F4
# ITS1F CTTGGTCATTTAGAGGAAGTAA
# ITS4 TCCTCCGCTTATTGATATGC
# ITS1F(revcomp) TTACTTCCTCTAAATGACCAAG
# ITS4(revcomp) GCATATCAATAAGCGGAGGA
PRIMERS_F[3]="CTTGGTCATTTAGAGGAAGTAA"
PRIMERS_FRC[3]="TTACTTCCTCTAAATGACCAAG"
PRIMERS_R[3]="TCCTCCGCTTATTGATATGC"
PRIMERS_RRC[3]="GCATATCAATAAGCGGAGGA"
PRODUCT_LENGTH_MIN[3]="400"
PRODUCT_LENGTH_MAX[3]="900"
# ITS54
# ITS5 GGAAGTAAAAGTCGTAACAAGG
# ITS5(revcomp) CCTTGTTACGACTTTTACTTCC
PRIMERS_F[4]="GGAAGTAAAAGTCGTAACAAGG"
PRIMERS_FRC[4]="CCTTGTTACGACTTTTACTTCC"
PRIMERS_R[4]="TCCTCCGCTTATTGATATGC"
PRIMERS_RRC[4]="GCATATCAATAAGCGGAGGA"
PRODUCT_LENGTH_MIN[4]="400"
PRODUCT_LENGTH_MAX[4]="900"
# ITS14
# ITS1 TCCGTAGGTGAACCTGCGG
# ITS1(revcomp) CCGCAGGTTCACCTACGGA
PRIMERS_F[5]="TCCGTAGGTGAACCTGCGG"
PRIMERS_FRC[5]="CCGCAGGTTCACCTACGGA"
PRIMERS_R[5]="TCCTCCGCTTATTGATATGC"
PRIMERS_RRC[5]="GCATATCAATAAGCGGAGGA"
PRODUCT_LENGTH_MIN[5]="400"
PRODUCT_LENGTH_MAX[5]="900"
# LepCOI
# LepF1 ATTCAACCAATCATAAAGATAT
# LepF1 (revcomp) ATATCTTTATGATTGGTTGAAT
# LepR1 TAAACTTCTGGATGTCCAAAAA
# LepR1 (revcomp) TTTTTGGACATCCAGAAGTTTA
PRIMERS_F[6]="ATTCAACCAATCATAAAGATAT"
PRIMERS_FRC[6]="ATATCTTTATGATTGGTTGAAT"
PRIMERS_R[6]="TAAACTTCTGGATGTCCAAAAA"
PRIMERS_RRC[6]="TTTTTGGACATCCAGAAGTTTA"
PRODUCT_LENGTH_MIN[6]="400"
PRODUCT_LENGTH_MAX[6]="1000"
# LSU
# LR0R ACCCGCTGAACTTAAGC
# LR0R (revcomp) GCTTAAGTTCAGCGGGT
# LR5 TCCTGAGGGAAACTTCG
# LR5 (revcomp) CGAAGTTTCCCTCAGGA
PRIMERS_F[7]="ACCCGCTGAACTTAAGC"
PRIMERS_FRC[7]="GCTTAAGTTCAGCGGGT"
PRIMERS_R[7]="TCCTGAGGGAAACTTCG"
PRIMERS_RRC[7]="CGAAGTTTCCCTCAGGA"
PRODUCT_LENGTH_MIN[7]="800"
PRODUCT_LENGTH_MAX[7]="1400"
# TEF
# 728 CATCGAGAAGTTCGAGAAGG
# 728 (revcomp) CCTTCTCGAACTTCTCGATG
# 1567 ACHGTRCCRATACCACCRATCTT
# 1567 (revcomp) AAGATYGGTGGTATYGGYACDGT
PRIMERS_F[8]="CATCGAGAAGTTCGAGAAGG"
PRIMERS_FRC[8]="CCTTCTCGAACTTCTCGATG"
PRIMERS_R[8]="ACHGTRCCRATACCACCRATCTT"
PRIMERS_RRC[8]="AAGATYGGTGGTATYGGYACDGT"
PRODUCT_LENGTH_MIN[8]="550"
PRODUCT_LENGTH_MAX[8]="1000"
length=9 # not doing 16S yet
readarray -t barcodes < barcodes.txt
readarray -t primers < primers.txt
i=0
#for folder in barcode*
for folder in "${barcodes[@]}"
do
	if [ -d "$folder" ]
	then
		cd $folder
		# a) join files together into one .fastq
		#rm all*.fastq
		cat *.fastq > all.fastq
		# b) generate quality plots
		NanoPlot --fastq all.fastq --plots hex
		# c) quality filtering and removing barcode adapters
		# Q10, minimum length 400 bp, trim leading 24 bases (adapter)
		NanoFilt -q 10 -l 200 --maxlength 1400 --headcrop 50 all.fastq > all_filt.fastq
		# d) convert to .fasta format
		sed -n '1~4s/^@/>/p;2~4p' all_filt.fastq > all_filt.fasta
		# e) shorten names of sequences to just 8 characters
		sed 's/.//10g; n' all_filt.fasta > all_filt_concatenated.fasta
		# f) reorient and trim primer sequences (allows for 10% error rate)
		j=""
		for k in "${!ARRAY_FOLDERS[@]}"; do
			if [[ "${ARRAY_FOLDERS[k]}" == "${primers[i]}" ]]; then
				j=$k;
			fi
		done
		if [[ "$j" == "" ]] || [[ "$j" -ge "$length" ]]; then
			echo "no primer named ${primers[i]}"
		else
			mkdir ${ARRAY_FOLDERS[j]}
			cd ${ARRAY_FOLDERS[j]}
			cp ../../input .
			cp ../../metabarcoding_consensus_files.R .
			cp ../../cons.denovo.py .
			cp ../../consensus_seq_and_score.denovo.R .
			cp ../../ReverseComplement.pl .
			cp ../../dnadist .
			# forward F reverse RRC
			cutadapt -g ${PRIMERS_F[j]} -a ${PRIMERS_RRC[j]} -o forward_seqs.fasta -m ${PRODUCT_LENGTH_MIN[j]} -M ${PRODUCT_LENGTH_MAX[j]} --overlap 15 ../all_filt_concatenated.fasta --untrimmed-output=all_filt_noF.fasta
			# forward R reverse FRC
			cutadapt -g ${PRIMERS_R[j]} -a ${PRIMERS_FRC[j]} -o revcomp_seqs.fasta -m ${PRODUCT_LENGTH_MIN[j]} -M ${PRODUCT_LENGTH_MAX[j]} --overlap 15 all_filt_noF.fasta --discard-untrimmed
			perl ReverseComplement.pl revcomp_seqs.fasta revcomp_seqs_revcomp.fasta
			cat forward_seqs.fasta revcomp_seqs_revcomp.fasta > all_filt_reorient.fasta
			myfilesize1=$(wc -l "all_filt_reorient.fasta" | awk '{print $1}')
			if (( $myfilesize1 <6))
			then
				cd ..
				rm -R ${ARRAY_FOLDERS[j]}
			else
#				/usr/local/bin/mafft --mpi --globalpair --thread 64 --op 0.5 --gop -0.5 --phylipout --quiet all_filt_reorient.fasta > all_filt_reorient.maf				
				/usr/local/bin/mafft --mpi --globalpair --thread 64 --threadtb 20 --threadit 20 --op 0.5 --gop -0.5 --phylipout --quiet all_filt_reorient.fasta > all_filt_reorient.maf
				myfilesizemafft=$(wc -l "all_filt_reorient.maf" | awk '{print $1}')
				if (( $myfilesizemafft == 0 ))
				then
					echo "global alignment for ${ARRAY_FOLDERS[j]} failed"
				else
					# phylip step
					# run the command file
					./dnadist < input
					mv outfile all_filt_reorient.dist
					echo "cluster(phylip=./all_filt_reorient.dist, method=average, cutoff=.2)" | ~/mothur/mothur
					echo "bin.seqs(list=./all_filt_reorient.an.list, fasta=./all_filt_reorient.fasta)" | ~/mothur/mothur
					# then get consensus sequences for each otu
					Rscript ./metabarcoding_consensus_files.R
					otus=$(ls Otu*.fasta | wc | awk '{print $1}')
					if (( $otus == 0 ))
					then
						cd ..
						rm -R ${ARRAY_FOLDERS[j]}
					else
						for otu_file in Otu*.fasta
						do
							# changed the following lines
							# default retree is 2
							head -n 1000 $otu_file > ${otu_file%.*}.concat.fasta
							# changes made 12/20/22
#							#/usr/local/bin/ginsi --clustalout --inputorder --thread 64 --quiet ${otu_file%.*}.concat.fasta > all_filt_denovo.maf
							/usr/local/bin/ginsi --clustalout --inputorder --thread 64 --threadtb 20 --threadit 20 --quiet ${otu_file%.*}.concat.fasta > all_filt_denovo.maf
							python cons.denovo.py
							mv all_filt_denovo.maf $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf
							Rscript consensus_seq_and_score.denovo.R $otu_file ${otu_file%.*}	
							mv all_filt_denovo.maf.karuna.cons $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.karuna.cons
							mv all_filt_denovo.maf.cons.data $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.cons.data
							mv all_filt_denovo.maf.karuna.cons.ungap $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.karuna.cons.ungap
							mv ${ARRAY_FOLDERS[j]}.cons.denovo.python.fasta $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.cons.denovo.python.fasta
							mv ${ARRAY_FOLDERS[j]}.cons.denovo.python.data $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.cons.denovo.python.data
						done
						rm ./input
						rm ./metabarcoding_consensus_files.R
						rm ./cons.denovo.py
						rm ./consensus_seq_and_score.denovo.R
						rm ./ReverseComplement.pl
						rm ./dnadist
						cd ..
					fi
				fi
			fi
			cd ..
		fi
	else
		echo "folder $folder does not exist"
	fi
	i=$((i+1))
done
# how do we want this organized ? for each barcode,
# a folder of results with the consensus sequences
# of all the otus in a separate folder for each type of 
mkdir results
cp ./barcode*/*/barcode*.*.Otu*.cons.denovo.python.fasta results
cd results
# added another file
cat *.fasta > all_seqs.fasta
cat $(echo $(ls | grep -P "barcode*.*.Otu[0]*[12345][[\.]]*")) > top_seqs.fasta
### blast
blastn -query top_seqs.fasta -db ~/blastdb/nt -max_target_seqs 1 -num_threads 64 -outfmt "6 qseqid sacc pident length mismatch evalue bitscore staxids sscinames stitle" -out blast.out
