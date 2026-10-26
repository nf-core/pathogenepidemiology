# Reference module replacement tests

From the repository root, verify the pinned launcher and run nf-test:

```bash
pixi run nextflow -version
pixi run nf-test version
pixi run nf-test test \
    -c subworkflows/local/prepare_references/tests/nf-test.config \
    modules/nf-core/wget/tests/main.nf.test \
    modules/nf-core/gunzip/tests/main.nf.test \
    subworkflows/local/prepare_references/tests/main.nf.test --ci
```

Both version commands must report Nextflow 26.10.3. Docker is required. The
focused runner loads the upstream WGET test dependency `nft-utils@0.0.9`
without changing the pipeline-wide runner or the installed canonical tests.
Canonical module tests fetch only their small nf-core test fixtures.

The subworkflow tests use `-stub` and synthetic, non-resolving URLs: no reference
genomes are downloaded. They check the WGET tuple interface, query-only GUNZIP,
query metadata and FAI, compressed host path output, publication under
`<outdir>/download`, and rejection of colliding URL basenames before tasks run.
They test wiring/publication, not FASTA content; the canonical tests exercise
the real WGET and GUNZIP commands as well as their stubs.
