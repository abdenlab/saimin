#!/usr/bin/env python

import pathlib
import sys


def pysam_run(file: pathlib.Path, chrom: str, start: int, stop: int):
    import pysam
    import pandas as pd

    ref_names = []
    starts = []
    ends = []
    names = []
    cigars = []
    seqs = []
    quals = []

    bam = pysam.AlignmentFile(str(file), "rb")

    for read in bam.fetch(chrom, start, stop):
        ref_names.append(read.reference_name)
        starts.append(read.reference_start)
        ends.append(read.reference_end)
        names.append(read.query_name)
        cigars.append(read.cigarstring)
        seqs.append(read.query_sequence)
        quals.append("".join(chr(ch + 33) for ch in read.query_qualities))

    return pd.DataFrame(
        {
            "ref_names": ref_names,
            "starts": starts,
            "ends": ends,
            "names": names,
            "cigars": cigars,
            "seqs": seqs,
            "quals": quals,
        }
    ).astype(
        {
            "ref_names": "category",
        }
    )


def saimin_ipc(file: pathlib.Path, chr: str, start: int, stop: int):
    import saimin

    with saimin.BamReader(str(file)) as reader:
        ipc = reader.fetch(chr, start, stop)

    return ipc


def saimin_pandas(file: pathlib.Path, chr: str, start: int, stop: int):
    import pyarrow
    from io import BytesIO

    ipc = saimin_ipc(file, chr, start, stop)
    df = pyarrow.ipc.open_file(BytesIO(ipc)).read_pandas()
    return df


def saimin_polars(file: pathlib.Path, chr: str, start: int, stop: int):
    import polars

    ipc = saimin_ipc(file, chr, start, stop)
    df = polars.read_ipc(ipc)
    return df


def main():
    file = pathlib.Path(__file__).parent / "data" / "example.bam"
    region = ("chr2", 1, 27_000_000)

    cmd = {
        "pysam": pysam_run,
        "saimin_pandas": saimin_pandas,
        "saimin_polars": saimin_polars,
        "saimin_ipc": saimin_ipc,
    }[sys.argv[1]]

    df = cmd(file, *region)
    print(df)


if __name__ == "__main__":
    main()
