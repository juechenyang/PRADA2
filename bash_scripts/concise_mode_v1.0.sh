echo 'run STAR two-pass way mapping ********************************************************************'

#adjust the file size limit for writing
ulimit -n 4096


#STAR command to perform alignment
#STAR 	--runMode genomeGenerate \
#     	--genomeDir ./genomeIndices \
#     	--genomeFastaFiles ./reference_data/reference.fasta \
#     	--sjdbOverhang 100 \
#     	--sjdbGTFfile ./reference_data/annotation.RNA-SeQC.gtf \
#     	--runThreadN $6 \
#     	--genomeChrBinNbits 16
#run STAR two-pass way mapping
if STAR 	--genomeDir ${10} \
	--readFilesIn $1 $2 \
	--chimOutJunctionFormat 1 \
	--outReadsUnmapped Fastx \
	--runThreadN $6 \
	--chimSegmentMin 12 \
	--chimJunctionOverhangMin 12 \
	--chimSegmentReadGapMax 3 \
	--alignSJDBoverhangMin 10 \
	--alignMatesGapMax 100000 \
	--alignIntronMax 100000 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--twopassMode Basic \
	--outFileNamePrefix $3/bam_results/ \
	--genomeLoad NoSharedMemory \
	--sjdbOverhang 100 \
	--outSAMstrandField intronMotif \
	--outSAMattributes NH HI NM MD AS XS \
	--outSAMunmapped Within \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 50000000000;
then
echo 'create index for bam file ********************************************************************'
sorted_bam=$3/bam_results/Aligned.sortedByCoord.out.bam
samtools index $sorted_bam



if [ "$9" = "True" ]
then
    #get variable names
    filename=$3/bam_results/Aligned.toTranscriptome.out.bam
    actual_name=${filename:0:${#filename}-4}
    marked_bam=$actual_name.marked_dup.bam



    # mark duplicate read of bam
    echo 'processing mark duplicates ********************************************************************'
    java -jar ${10}/tools/picard.jar MarkDuplicates \
          I=$filename \
          O=$marked_bam \
          M=$3/inter_results/marked_dup_metrics.txt \
          VALIDATION_STRINGENCY=LENIENT
fi



# whether to keep original bam
if [ "$7" = "False" ]
then
    rm $3/bam_results/Aligned.toTranscriptome.out.bam
fi

# whether to keep intermediate file
if [ "$5" = "False" ]
then
    rm -rf $3/inter_results
fi

# whether to perform RNA-SeQC
if [ "$4" = "True" ]
then
    # build index of bam
    echo 'building index for aligned BAM file ********************************************************************'
    samtools index $marked_bam



    #switch to Java 1.7
    export PATH="/home/CBBI/yangj8/Java/1.7/jdk1.7.0_80/bin:$PATH"

    echo 'performing quality control ********************************************************************'
    # generate QC report
    java -Xmx8g \
         -jar ${10}/tools/RNA-SeQC_v1.1.8.jar \
         -ttype 2 \
         -r ${10}/reference_data/reference.fasta \
         -o $3/QC/ \
         -s "$8|$marked_bam|Disc" \
         -t ${10}/reference_data/annotation.RNA-SeQC.gtf
fi
else
exit 1
    
fi



