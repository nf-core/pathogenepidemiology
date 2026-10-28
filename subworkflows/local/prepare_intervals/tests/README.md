# PREPARE_INTERVALS tests

Run from the repository root, using the Pixi-pinned Nextflow launcher:

```bash
NFT_WORKDIR=/tmp/opencode/nft-intervals pixi run nf-test test \
  subworkflows/local/prepare_references/tests/main.nf.test \
  subworkflows/local/prepare_intervals/tests/main.nf.test \
  subworkflows/local/gatk_MOI/tests/concat.nf.test \
  subworkflows/local/gatk_MOI/tests/main.nf.test \
  -c subworkflows/local/prepare_intervals/tests/nf-test.config
```

These tests use the tracked `assets/Pf3D7_core.bed` (including reordered and
single-contig subsets), shell fakes and stubbed analysis processes to
check BED validation and reuse, scatter wiring, VQSR paths and concat ordering.
Reference setup additionally checks value-channel types and two-sample reuse
with a host present or absent; download/decompression commands are shell fakes.
They do not run real GATK analysis, verify coordinate correctness, test HTTP
transfers or check real command handling of paths with spaces. Logs and
artifacts remain in `NFT_WORKDIR`.

[GATK 4.7.0.0 `SplitIntervals.onTraversalStart`](https://github.com/broadinstitute/gatk/blob/4.7.0.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/SplitIntervals.java) scatters once with
`--scatter-count 1`, then groups that list by contig with
`--dont-mix-contigs true`, independent of BED input order. The global duplicate
check is defensive; its failure test copies a fake output before conversion.
Concat retains lexical filename order; no natural-order sort was added.

Resume caching is not asserted by this suite: nf-test recreates each test's
Nextflow cache when rerunning it, so adding `-resume` alone does not verify reuse.
