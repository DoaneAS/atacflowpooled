#!/bin/bash
p1=$1

if [ -z "$p1" ]
then
    echo "This removes sup chromosomes from BAM files"
    echo "$(basename $0) <bamFile>"
    echo "First input is the bamFile"
    exit
fi

MYSLOTS=$2

if [ -z ${NSLOTS+x}  ]
then
    NSLOTS=$MYSLOTS
fi

samtools index ${p1}

o1a=$(echo $(basename ${p1}) | sed -r 's/\.bam$/.nosup.bam/g')

samtools idxstats ${p1} | cut -f 1 |  head -n 23 | xargs samtools view -@ ${NSLOTS} -b ${p1} > ${o1a}
