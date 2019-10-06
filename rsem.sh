if rsem-calculate-expression -p 8 --paired-end \
    --bam \
    --estimate-rspd \
    --append-names \
    --output-genome-bam \
    $1 \
    $2 $3
then
exit 0
else
exit 1
fi

