#!/usr/bin/env bash

#adjust the file size limit for writing
ulimit -n 4096


#STAR command to perform alignment
STAR 	--runMode genomeGenerate \
     	--genomeDir $1 \
     	--genomeFastaFiles $2 \
     	--sjdbOverhang 100 \
     	--sjdbGTFfile $3 \
     	--runThreadN $4 \
     	--genomeChrBinNbits 16