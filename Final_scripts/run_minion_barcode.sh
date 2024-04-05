#!/bin/bash
rm ./metabarcode_output.log
sh ./Master_metabarcoding_script_mpi.sh >> ./metabarcode_output.log 2>&1
