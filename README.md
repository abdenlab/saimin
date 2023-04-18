# bram

proof of concept. read short read sequence data stored in sam/bam files as arrow ipc (without htslib).

# usage

```python
import polars
import bram

with bram.BamReader("example.bam") as reader:
    ipc = reader.fetch("chr2", 0, 100_000)
df = polars.read_ipc(ipc)
df

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
