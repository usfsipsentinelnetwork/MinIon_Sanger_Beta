MinIon_Sanger_Beta - for publication of sequences from Gnathotrupes

(Very) basic how-to-use (better annotated versions of scripts/ with flags/options will be available in MIMP package)

0. Install dependencies (see header of Master_metabarcoding_script_mpi.sh)

Linux
# EMBL_EBI
# MAFFT
# emboss embassy package phylip implementation (ednadist)
# python
# Databases
# mothur
# NanoFilt
# cutadapt
# R

1. Edit barcodes.txt
2. Edit primers.txt
3. Change any parameters (ie cutoff for cutadapt, nanofilt, etc.)

4. run Master_metabarcoding_script_mpi.sh (you can use run_minion_barcode.sh)
