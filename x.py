#!/usr/bin/env python

import sys

import polars

import bram


def main():
    with bram.BamReader(sys.argv[1]) as reader:
        ipc = reader.fetch("chr1", 0, 100_000)
    df = polars.read_ipc(ipc)
    print(df.head())


if __name__ == "__main__":
    main()
