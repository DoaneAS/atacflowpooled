#!/bin/bash -l


bam=$1
bed=$2
Sample=$3


MYSLOTS=$4
if [ -z ${NSLOTS+x}  ]
then
    NSLOTS=$MYSLOTS
fi




###spack load samtools
#samtools index ${bam}
##gunzip ${Sample}.tn5.broadPeak.gz

##bedtools slop -i ${Sample}.tn5.broadPeak -g ${chrsz} -b 1000 > ${Sample}.slop1k.bed

nucleoatac run --bed ${bed} \
           --bam ${bam} \
           --fasta /athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
           --out ${Sample} \
           --write_all \
           --cores ${NSLOTS}




#nucleoatac run --bed "${TMPDIR}/${Sample}/${Sample}.slop1k.bed" --bam "${TMPDIR}/${Sample}/${Sample}.sorted.nodup.noM.black.bam" \
#           --fasta $TMPDIR/$RGEN --out  $TMPDIR/${Sample}/nucATAC/${Sample} \
#           --write_all --cores ${NSLOTS}
