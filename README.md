# saimin 🍜

proof of concept. read short read sequence data stored in sam/bam files as arrow ipc (without htslib).

## usage

```python
import polars
import saimin

with saimin.BamReader("./example.bam") as reader:
    ipc = reader.fetch("chr2", 1, 1_000_000)

df = polars.read_ipc(ipc)
print(df)

# shape: (655, 7)
# ┌──────┬────────┬────────┬───────────────────────────────────┬───────┬───────────────────────────────────┬───────────────────────────────────┐
# │ ref  ┆ start  ┆ end    ┆ name                              ┆ cigar ┆ seq                               ┆ qual                              │
# │ ---  ┆ ---    ┆ ---    ┆ ---                               ┆ ---   ┆ ---                               ┆ ---                               │
# │ cat  ┆ i32    ┆ i32    ┆ str                               ┆ str   ┆ str                               ┆ str                               │
# ╞══════╪════════╪════════╪═══════════════════════════════════╪═══════╪═══════════════════════════════════╪═══════════════════════════════════╡
# │ chr2 ┆ 10011  ┆ 10047  ┆ SOLEXA-1GA-2_2_FC20EMB:5:5:440:4… ┆ 36M   ┆ CACCAGACCCACACACCAAACCCACACACACA… ┆ G3<@3<'-**)/&4",,'%'**%,&)(%#'#…  │
# │ chr2 ┆ 10011  ┆ 10047  ┆ SOLEXA-1GA-2_2_FC20EMB:5:116:68:… ┆ 36M   ┆ CACCACACCCACAGACCACACCCACACACACA… ┆ "  ,$'")$0%!'% ,/++4,->+?16'…     │
# │ chr2 ┆ 10018  ┆ 10054  ┆ SOLEXA-1GA-2_2_FC20EMB:5:296:593… ┆ 36M   ┆ CCCACACACCACACCCAGACACACACACACCC… ┆ 33GGGG27G423/1)G?&&14GGG9#%,!"'"… │
# │ chr2 ┆ 10023  ┆ 10059  ┆ SOLEXA-1GA-2_2_FC20EMB:5:11:343:… ┆ 36M   ┆ ACCCCACACCCACACACACCCACACACACACC… ┆ /-8GG:G,C&*(4(@+? "2)=*A55"7*4…   │
# │ …    ┆ …      ┆ …      ┆ …                                 ┆ …     ┆ …                                 ┆ …                                 │
# │ chr2 ┆ 989952 ┆ 989988 ┆ SOLEXA-1GA-2_2_FC20EMB:5:294:373… ┆ 36M   ┆ CACCTGTGAATTGTCCTGGAAATCCTGTTGAA… ┆ 9 1/$#?5G$26G$E74F57G+GG:9GGGGG…  │
# │ chr2 ┆ 992481 ┆ 992517 ┆ SOLEXA-1GA-2_2_FC20EMB:5:68:250:… ┆ 36M   ┆ TTCTGTAGCAGTGAATGAACAAAAGGAGCAAA… ┆ GGGGGG:F;5G/G40G50(3%,*-,$#&+$$!… │
# │ chr2 ┆ 998004 ┆ 998040 ┆ SOLEXA-1GA-2_2_FC20EMB:5:297:621… ┆ 36M   ┆ ATTATAGAATGGGGGAAAATATTTGCAAACTG… ┆ %-$E*G:GG,=GGGGGGGGGGGGGGGGGGGGG… │
# │ chr2 ┆ 998984 ┆ 999020 ┆ SOLEXA-1GA-2_2_FC20EMB:5:132:376… ┆ 36M   ┆ GAGGGTGACAGTGTATAGTCAAGAGTGTGGCC… ┆ "%/ $1A,!03)'/+3?5*<0:G>DGGCG;G2… │
# └──────┴────────┴────────┴───────────────────────────────────┴───────┴───────────────────────────────────┴───────────────────────────────────┘
```

## development

this project uses `maturin` and `hatch` for development, which can be installed with `pipx`.

```sh
# create a virtual env
hatch shell

# compile a development version of the package
maturin develop --release

# run the benchmark
hatch run bench

# Benchmark 1: ./bench.py pysam
#   Time (mean ± σ):     669.5 ms ±   8.8 ms    [User: 1200.9 ms, System: 363.7 ms]
#   Range (min … max):   660.8 ms … 690.3 ms    10 runs
#
# Benchmark 2: ./bench.py saimin_polars
#   Time (mean ± σ):     197.9 ms ±   1.7 ms    [User: 141.5 ms, System: 40.5 ms]
#   Range (min … max):   194.0 ms … 200.2 ms    14 runs
#
# Benchmark 3: ./bench.py saimin_pandas
#   Time (mean ± σ):     601.5 ms ±   5.5 ms    [User: 1153.5 ms, System: 352.1 ms]
#   Range (min … max):   597.6 ms … 616.0 ms    10 runs
#
# Benchmark 4: ./bench.py saimin_ipc
#   Time (mean ± σ):     122.0 ms ±   1.7 ms    [User: 90.5 ms, System: 19.9 ms]
#   Range (min … max):   119.0 ms … 124.8 ms    22 runs
#
# Summary
#   './bench.py saimin_ipc' ran
#     1.62 ± 0.03 times faster than './bench.py saimin_polars'
#     4.93 ± 0.08 times faster than './bench.py saimin_pandas'
#     5.49 ± 0.10 times faster than './bench.py pysam'
```
