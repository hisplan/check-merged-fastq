import os
import glob
import gzip
import re
import argparse
import pandas as pd
import ray


def parse_fastq_name(sample_name: str, read_type: str, fastqs: list) -> list:

    if not read_type in ["R1", "R2"]:
        raise Exception("read_type must be either R1 or R2.")

    results = []

    for fn in fastqs:
        pattern = f"{sample_name}_L(\d+?)_{read_type}_.*\.fastq\.gz"
        match = re.search(pattern, fn)
        if not match:
            raise Exception("Invalid FASTQ filename!")

        lane_name = match.group(1)

        results.append((read_type, lane_name, fn))

    return results


def parse_r1(sequence: str, ten_x_kit: str):
    if ten_x_kit == "v3":
        umi_len = 12
    elif ten_x_kit == "v2":
        umi_len = 10
    else:
        raise Exception("ten_x_kit must be either v2 or v3.")
    cb, umi = sequence[:16], sequence[16 : 16 + umi_len]
    poly_t = sequence[16 + umi_len :]
    return cb, umi, poly_t


def read_fastq_line1(filename: str):
    read_id = None
    sequence = None
    with gzip.open(filename, "rt") as fin:
        for line in fin:
            line = line.strip()
            if not read_id:
                # @A00333:373:HF27HDSX2:3:1101:1597:1000 1:N:0:NCCGTTCT+NCAATCCGTC
                read_id = line.split(" ")[0]
            elif not sequence:
                sequence = line
            else:
                break

    return read_id, sequence


def get_merged_read_id(cb: str, umi: str, poly_t: str, read_id: str) -> str:
    return f"@:{cb}:{umi}:{poly_t};{read_id[1:]}"


@ray.remote
def check(
    merged_fastq: str, cb: str, umi: str, poly_t: str, read_id: str, sequence: str
):
    # read_id: expected
    # sequence: expected

    # @:GTGTAACTCATACGAC:TCATATCAATGT:T;A00333:373:HF27HDSX2:2:1101:2899:1000 2:N:0:NCCGTTCT+NCAATCCGTC
    # CCGCTGCACAGGCTGCCTTCCAGAAGGTGGTGGCTGGAGTGGCCGCTGCCCTGGCTCACAAGTACCACTAAACCCCCTTTCCTGCTCTTGC
    # +
    # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFF

    # we're going to find by read_id e.g. A00333:373:HF27HDSX2:2:1101:2899:1000
    # remove `@` chracter
    query = read_id[1:]

    read_id_expected = get_merged_read_id(cb, umi, poly_t, read_id)
    sequence_expected = sequence

    # print(f"Searching for `{query}` in {merged_fastq}...")

    read_id_actual = None
    sequence_actual = None
    count = 0

    with gzip.open(merged_fastq, "rt") as fin:
        for line in fin:
            line = line.strip()
            if not read_id_actual:
                # @:GTGTAACTCATACGAC:TCATATCAATGT:T;A00333:373:HF27HDSX2:2:1101:2899:1000 2:N:0:NCCGTTCT+NCAATCCGTC
                read_id_actual = line.split(" ")[0]
            elif not sequence_actual:
                sequence_actual = line
                if query in read_id_actual:
                    if read_id_expected != read_id_actual:
                        # found, incorrect
                        return -1, merged_fastq, read_id_actual, sequence_actual

                    # found, looks okay
                    return 0, merged_fastq, read_id_actual, sequence_actual

            count += 1
            if count == 4:
                read_id_actual = None
                sequence_actual = None
                count = 0

    # not found
    return 1, None, None, None


def main(
    sample_name: str,
    dir_barcodes: str,
    dir_genomic: str,
    merged_fastqs: list,
    ten_x_kit: str,
):

    barcode_fastqs = glob.glob(os.path.join(dir_barcodes, "*.fastq.gz"))
    genomic_fastqs = glob.glob(os.path.join(dir_genomic, "*.fastq.gz"))

    if len(barcode_fastqs) == 0:
        raise Exception(f"No *.fastq.gz found in {dir_barcodes}")

    if len(genomic_fastqs) == 0:
        raise Exception(f"No *.fastq.gz found in {dir_genomic}")

    r1_fastqs = parse_fastq_name(
        sample_name=sample_name, read_type="R1", fastqs=barcode_fastqs
    )
    r2_fastqs = parse_fastq_name(
        sample_name=sample_name, read_type="R2", fastqs=genomic_fastqs
    )

    column_names = ["read_type", "lane_num", "fastq"]
    df_r1 = pd.DataFrame(r1_fastqs, columns=column_names)
    df_r2 = pd.DataFrame(r2_fastqs, columns=column_names)

    print(df_r1.to_markdown(tablefmt="github"))
    print()
    print(df_r2.to_markdown(tablefmt="github"))

    # check equal lane number
    lane_numbers = df_r1.lane_num.unique()

    for lane_num in lane_numbers:

        print()
        print("Lane:", lane_num)

        r1_fastq = df_r1[df_r1.lane_num == lane_num].fastq.values[0]
        r2_fastq = df_r2[df_r2.lane_num == lane_num].fastq.values[0]

        print("R1:", r1_fastq)
        print("R2:", r2_fastq)

        read_id_r1, sequence_r1 = read_fastq_line1(r1_fastq)
        cb, umi, poly_t = parse_r1(sequence=sequence_r1, ten_x_kit=ten_x_kit)

        read_id_r2, sequence_r2 = read_fastq_line1(r2_fastq)

        if read_id_r1 != read_id_r2:
            print(read_id_r1)
            print(read_id_r2)
            raise Exception("read_id from R1 and R2 do not match!")

        futures = []
        for merged_fastq in merged_fastqs:
            print(f"Searching for `{read_id_r1[1:]}` in {merged_fastq}...")
            future = check.remote(
                merged_fastq=merged_fastq,
                cb=cb,
                umi=umi,
                poly_t=poly_t,
                read_id=read_id_r1,  # read_id_r2 is the same as read_id_r1
                sequence=sequence_r2,
            )
            futures.append(future)

        keep_check = True
        while keep_check:
            ready, not_ready = ray.wait(futures, num_returns=len(futures), timeout=5)

            print("Ready:", len(ready), "Not Ready:", len(not_ready))

            # see if any of the ready one returns 0 or -1 (i.e. found)
            for future in ready:
                return_code, filename, read_id_actual, sequence_actual = ray.get(future)
                # found, okay
                if return_code == 0:
                    print(f"Found in `{filename}`")
                    print("Looks okay!")
                    keep_check = False
                # found, not okay
                elif return_code != 1:
                    print(f"Found in `{filename}`")
                    print("> Expected")
                    print(get_merged_read_id(cb, umi, poly_t, read_id_r1))
                    print(sequence_r2)
                    print()
                    print("> Actual")
                    print(read_id_actual)
                    print(sequence_actual)

                    ray.shutdown()
                    exit(1)

            # exit if all ready
            if len(not_ready) == 0:
                keep_check = False


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample",
        action="store",
        dest="sample_name",
        help="[SAMPLE NAME]_S1_L00[LANE NUMBER]_[READ TYPE]_001.fastq.gz",
        required=True,
    )

    parser.add_argument(
        "--barcode",
        action="store",
        dest="dir_barcode",
        help="Directory containing barcode FASTQ files",
        required=True,
    )

    parser.add_argument(
        "--genomic",
        action="store",
        dest="dir_genomic",
        help="Directory containing genomic FASTQ files",
        required=True,
    )

    parser.add_argument(
        "--kit", action="store", dest="ten_x_kit", help="either v2 or v3", required=True
    )

    parser.add_argument(
        "--chunk-prefix",
        action="store",
        dest="chunk_prefix",
        default="chunk",
        help="Prefix used when splitting your merged FASTQ file into multiple chunks. Include the directory name.",
    )

    parser.add_argument(
        "--threads",
        action="store",
        type=int,
        dest="num_threads",
        default=20,
        help="Number of parallel tasks to be used to speed up the search",
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    ray.init(num_cpus=params.num_threads)

    # sample_name = "2653_blood2_CX3CR1_CCR2_IGO_12104_39_S13"
    # ten_x_kit = "v3"
    # dir_barcodes = "/fscratch/chunj/seqc-sort/barcode"
    # dir_genomic = "/fscratch/chunj/seqc-sort/genomic"

    # sample_name = "1890_RA19_10_13_DAPI_Low_IGO_10875_23"
    # ten_x_kit = "v3"
    # dir_barcodes = "/fast/chunj/nucseq/fastq/barcode"
    # dir_genomic = "/fast/chunj/nucseq/fastq/genomic"

    merged_fastqs = glob.glob(f"{params.chunk_prefix}-*.fastq.gz")

    if len(merged_fastqs) == 0:
        raise Exception("No chunked merged fastq file found.")

    main(
        sample_name=params.sample_name,
        dir_barcodes=params.dir_barcode,
        dir_genomic=params.dir_genomic,
        merged_fastqs=merged_fastqs,
        ten_x_kit=params.ten_x_kit,
    )
