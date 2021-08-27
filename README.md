# check-merged-fastq

Check whether a merged FASTQ file generated by SEQC is correct or not.

## Set Up

```
conda create -n check python=3.8 pip
conda activate check
pip install -r requirements.txt
```

## How to Run

Download barcode/genomic FASTQ files and merged FASTQ file from AWS S3:

```bash
$ cd workspace
$ ./download.sh
USAGE: download.sh [options]
    -s  S3 URI (e.g. s3://dp-lab-data/collaborators/pi/prj/blood2_CX3CR1_CCR2)
    -d  local destination (e.g. blood2_CX3CR1_CCR2)

$ ./download.sh \
  -s s3://dp-lab-data/collaborators/pi/prj/blood2_CX3CR1_CCR2 \
  -d blood2_CX3CR1_CCR2
```

```
workspace
├── blood2_CX3CR1_CCR2-HPC
│   ├── blood2_CX3CR1_CCR2_merged.fastq.gz
│   ├── barcode
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R1_001.fastq.gz
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R1_001.fastq.gz
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R1_001.fastq.gz
│   │   └── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R1_001.fastq.gz
│   └── genomic
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R2_001.fastq.gz
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R2_001.fastq.gz
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R2_001.fastq.gz
│       └── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R2_001.fastq.gz
```

NOTE: If you want to test the merged FASTQ file generated by running SEQC on HPC, make sure you overwrite the one downloaded from AWS with the one you generated from HPC.

Split the merged FASTQ files into multiple chunks in order to speed up the checking process:

```bash
$ cd ..
$ ./split.sh
USAGE: split.sh [options]
    -w  working directory (e.g. workspace/blood2_CX3CR1_CCR2)
    -m  merged FASTQ filename only, not the whole path (e.g. 2653_blood2_CX3CR1_CCR2_IGO_12104_39_merged.fastq.gz)
    -n  number of chunks to generate (default=20)

$ ./split.sh
  -w workspace/blood2_CX3CR1_CCR2 \
  -m 2653_blood2_CX3CR1_CCR2_IGO_12104_39_merged.fastq.gz
```

```
workspace
├── blood2_CX3CR1_CCR2-HPC
│   ├── blood2_CX3CR1_CCR2_merged.fastq.gz
│   ├── chunk-001.fastq.gz
│   ├── chunk-002.fastq.gz
│   ├── chunk-003.fastq.gz
│   ├── chunk-004.fastq.gz
│   ├── chunk-005.fastq.gz
│   ├── chunk-006.fastq.gz
│   ├── chunk-007.fastq.gz
│   ├── chunk-008.fastq.gz
│   ├── chunk-009.fastq.gz
│   ├── chunk-010.fastq.gz
│   ├── chunk-011.fastq.gz
│   ├── chunk-012.fastq.gz
│   ├── chunk-013.fastq.gz
│   ├── chunk-014.fastq.gz
│   ├── chunk-015.fastq.gz
│   ├── chunk-016.fastq.gz
│   ├── chunk-017.fastq.gz
│   ├── chunk-018.fastq.gz
│   ├── chunk-019.fastq.gz
│   ├── chunk-020.fastq.gz
│   ├── barcode
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R1_001.fastq.gz
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R1_001.fastq.gz
│   │   ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R1_001.fastq.gz
│   │   └── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R1_001.fastq.gz
│   └── genomic
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R2_001.fastq.gz
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R2_001.fastq.gz
│       ├── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R2_001.fastq.gz
│       └── 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R2_001.fastq.gz
```

Run the checking process:

```bash
$ python check_validity.py \
    --sample 2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13 \
    --barcode workspace/blood2_CX3CR1_CCR2/barcode \
    --genomic workspace/blood2_CX3CR1_CCR2/genomic \
    --chunk-prefix workspace/blood2_CX3CR1_CCR2/chunk \
    --kit v3 \
    --threads=20
```

NOTE:

- `--sample` should be set to the `SAMPLE NAME` prefix of the FASTQ file name: `[SAMPLE NAME]_S1_L00[LANE NUMBER]_[READ TYPE]_001.fastq.gz`
- For example, if your FASTQ file name is `2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R1_001.fastq.gz`, then `--sample` should be set to `2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13`.

Below is the example output of the search. If your merged FASTQ file is invalid, you will see "expected" vs. "actual", and the exit code will be 1. If your merged FASTQ file looks okay, you will see "looks okay", and the exit code will be 0.

```
|   | read_type | lane_num | fastq                                                                               |
|---|-----------|----------|-------------------------------------------------------------------------------------|
| 0 | R1        |      001 | workspace/.../barcode/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R1_001.fastq.gz |
| 1 | R1        |      002 | workspace/.../barcode/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R1_001.fastq.gz |
| 2 | R1        |      003 | workspace/.../barcode/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R1_001.fastq.gz |
| 3 | R1        |      004 | workspace/.../barcode/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R1_001.fastq.gz |

|   | read_type | lane_num | fastq                                                                               |
|---|-----------|----------|-------------------------------------------------------------------------------------|
| 0 | R2        |      001 | workspace/.../genomic/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R2_001.fastq.gz |
| 1 | R2        |      002 | workspace/.../genomic/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L002_R2_001.fastq.gz |
| 2 | R2        |      003 | workspace/.../genomic/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L003_R2_001.fastq.gz |
| 3 | R2        |      004 | workspace/.../genomic/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L004_R2_001.fastq.gz |

Lane: 001
R1: workspace/.../barcode/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R1_001.fastq.gz
R2: workspace/.../genomic/2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13_L001_R2_001.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-001.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-002.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-003.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-004.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-005.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-006.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-007.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-008.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-009.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-010.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-011.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-012.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-013.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-014.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-015.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-016.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-017.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-018.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-019.fastq.gz
Searching for `A00333:373:HF27HDSX2:1:1101:1127:1000` in workspace/.../chunk-020.fastq.gz
Ready: 0 Not Ready: 20
Ready: 0 Not Ready: 20
Ready: 1 Not Ready: 19
Found in `workspace/.../chunk-001.fastq.gz`
> Expected
@:GATGATCAGTAGCCAG:ATTCACTTTATG:G;A00333:373:HF27HDSX2:1:1101:1127:1000
GCATCAGCCTAGAGCAGGACAAGCCACGTCAGCCAGCTCTGATTTGACTGAGAAACTCTGCCTCAAAGAATAAGGCAGAGCAATCAAGGAT

> Actual
@:ACCAAACCACATGGTT:CATCTATAAGGT:G;A00333:373:HF27HDSX2:1:1101:1127:1000
GCATCAGCCTAGAGCAGGACAAGCCACGTCAGCCAGCTCTGATTTGACTGAGAAACTCTGCCTCAAAGAATAAGGCAGAGCAATCAAGGAT
```

```bash
echo $?
1
```
