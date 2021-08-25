#!/bin/bash -e

# no trailing slash

usage()
{
cat << EOF
USAGE: `basename $0` [options]
    -s  S3 URI (e.g. s3://dp-lab-data/collaborators/aboire/ImmunologicalDeterminantsLeptomeningealMet/blood2_CX3CR1_CCR2)
    -d  destination (e.g. blood2_CX3CR1_CCR2)
EOF
}

while getopts "s:d:h" OPTION
do
    case $OPTION in
        s) src=$OPTARG ;;
        d) dst=$OPTARG ;;
        h) usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$src" ] || [ -z "$dst" ]
then
    usage
    exit 1
fi

aws s3 sync --exclude="*" --include="*_merged.fastq.gz" ${src}/seqc-results/ ${dst}/
aws s3 sync ${src}/barcode/ ${dst}/barcode/
aws s3 sync ${src}/genomic/ ${dst}/genomic/
