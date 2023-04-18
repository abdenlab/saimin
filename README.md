# bram

proof of concept. read short read sequence data stored in sam/bam files as arrow ipc (without htslib).

# usage

```python
import polars
import bram

with bram.BamReader("example.bam") as reader:
    ipc = reader.fetch("chr1", 0, 100_000)
df = polars.read_ipc(ipc)
df

# shape: (5, 7)
# ┌────────┬───────┬───────┬────────────────────┬───────┬──────────────┬────────────────┐
# │ ref_id ┆ start ┆ end   ┆ name               ┆ cigar ┆ seq          ┆ qual           │
# │ ---    ┆ ---   ┆ ---   ┆ ---                ┆ ---   ┆ ---          ┆ ---            │
# │ i32    ┆ i32   ┆ i32   ┆ str                ┆ str   ┆ str          ┆ str            │
# ╞════════╪═══════╪═══════╪════════════════════╪═══════╪══════════════╪════════════════╡
# │ 0      ┆ 10144 ┆ 10180 ┆ SOLEXAB:5:251:979… ┆ 36M   ┆ AACCAACCCTA… ┆ GGGG'!9$'"!…   │
# │ 0      ┆ 10147 ┆ 10183 ┆ SOLEXAB:5:102:214… ┆ 36M   ┆ CCCTCCTAACC… ┆ GAEGG*&+)"'%%… │
# │ 0      ┆ 10148 ┆ 10184 ┆ SOLEXAB:5:195:284… ┆ 36M   ┆ CCAACTAACCT… ┆                │
# │ 0      ┆ 10149 ┆ 10185 ┆ SOLEXAB:5:35:583:… ┆ 36M   ┆ CTAATAACCTA… ┆ GG6>7)+"%&&) … │
# │ 0      ┆ 10151 ┆ 10187 ┆ SOLEXAB:5:248:130… ┆ 36M   ┆ AACCACCTAAC… ┆ GBG/G)+7. +&…  │
# └────────┴───────┴───────┴────────────────────┴───────┴──────────────┴────────────────┘
```

# development

this project uses `maturin` and `hatch` for development, which can be installed with `pipx`.

```sh
# create a virtual env
hatch shell

# compile a development version of the package
maturin develop

# run the example script
./x.py data/example.bam
```
