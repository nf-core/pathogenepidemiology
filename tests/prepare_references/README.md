# Reference preparation integration tests

Run from the repository root:

```bash
python -m unittest discover -s tests/prepare_references -v
```

These tests use `pixi run nextflow` and local `wget`, `gzip`, `samtools` and `gatk`.
They execute the real installed nf-core modules, not process stand-ins or
`-stub-run`. Containers and Conda are disabled only in this isolated harness.
Synthetic FASTAs are served over localhost; no biological datasets are fetched.
Logs, traces, published outputs and task artifacts remain in the printed
`/tmp/opencode/prepare-references-*` directory.

The harness bypasses pipeline initialization/plugins and downstream analyses.
It checks both reference inputs independently, local plain/gzip combinations,
URL fallbacks, invalid scalar/file inputs with parameter validation disabled,
query metadata/index joins, URL-sensitive resume caching, and the `download`
publication layout. They also check dictionary generation, reusable value-channel
types and the optional empty-host value introduced on dev. Host output
remains a plain path channel, the interface consumed by `PREPROCESS_READS` and
`BBDUK_FILTER`.

Local FASTA basenames must start with an ASCII letter, digit or underscore and
contain only ASCII letters, digits, dots, underscores and hyphens. The suite
checks early rejection of spaces, shell metacharacters and leading hyphens for
both plain and gzip query/host files, before any reference tasks launch.

Installed with nf-core/tools 4.1.0, with unmodified upstream files and revisions
recorded in `modules.json`:

- `wget`: `8dfb9868efd0e39f27ae0fed1ed204dc8f0b2a13`; input
  `[meta, url, suffix]`, output `outfile: [meta, file]`. The filename is
  `${meta.id}.${suffix}`; the URL is an explicit cache-key input.
- `gunzip`: `0902eac3012baaf4f9ab6513c8c55acc9353c96c`; input
  `[meta, archive]`, output `gunzip: [meta, file]`.

Upstream nf-test fixtures are included by the module installer and excluded
by the project's existing `modules/nf-core/**/tests/*` ignore rule. This Python
suite is an explicit focused integration check, not a new nf-test snapshot.
