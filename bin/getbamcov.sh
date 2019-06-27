#!/bin/bash

sbam=$1
outprefix=$2
sfactors=$3
MYSLOTS=$4

export OMP_NUM_THREADS=1

if [ -z ${NSLOTS+x}  ]
then
    NSLOTS=$MYSLOTS
fi


. /home/asd2007/.spackloads.sh

spack load samtools

samtools index ${sbam}

#source activate deepsge

bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
            --scaleFactor $sfactors \
            --maxFragmentLength 150 \
            -o ${outprefix}.deseqsf.bw \
            --centerReads \
            --extendReads \
            --numberOfProcessors ${NSLOTS}


##            --region chr3:186984410:188482691


##bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
##            --scaleFactor $sfactors \
##            --normalizeUsing CPM \
##            --maxFragmentLength 150 \
##            -o ${outprefix}.cpm.sizefactors.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}
##
##
##bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
##            --scaleFactor $sfactors \
##            --normalizeUsing BPM \
##            --maxFragmentLength 150 \
##            -o ${outprefix}.bpm.sizefactors.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}


#source deactivate
