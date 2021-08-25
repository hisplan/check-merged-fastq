#!/bin/bash

#merged_fastq="/fast/chunj/nucseq/seqc-results-pre-mrna/1890_RA19_10_13_DAPI_Low_IGO_10875_23_merged.fastq.gz"
# merged_fastq="blood2_CX3CR1_CCR2_merged.fastq.gz"
#merged_fastq="tmp/TEST_merged.fastq.gz"


# merged_fastq="2653_blood2_CX3CR1_CCR2_IGO_12104_39_merged.fastq.gz"
# workdir="workspace/blood2_CX3CR1_CCR2"

# default
n_chunks=20

usage()
{
cat << EOF
USAGE: `basename $0` [options]
    -w  working directory (e.g. workspace/blood2_CX3CR1_CCR2)
    -m  merged FASTQ filename only, not the whole path (e.g. 2653_blood2_CX3CR1_CCR2_IGO_12104_39_merged.fastq.gz)
    -n  number of chunks to generate (default=${n_chunks})
EOF
}

while getopts "w:m:n:h" OPTION
do
    case $OPTION in
        w) workdir=$OPTARG ;;
        m) merged_fastq=$OPTARG ;;
        n) n_chunks=$OPTARG ;;
        h) usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$workdir" ] || [ -z "$merged_fastq" ] || [ -z "$n_chunks" ]
then
    usage
    exit 1
fi

mkdir -p ${workdir}

outs=""
prefix="chunk"

for ((i = 1 ; i <= ${n_chunks} ; i++))
do
  chunk_num="`printf %03d $i`"
  outs+="-o ${workdir}/${prefix}-${chunk_num}.fastq.gz "
done

fastqsplitter -i ${workdir}/${merged_fastq} ${outs} -t ${n_chunks}
