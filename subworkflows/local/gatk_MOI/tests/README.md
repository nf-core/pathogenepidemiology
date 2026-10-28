# GATK_MOI lightweight tests

`main.nf.test` checks scatter wiring with analysis processes stubbed out.
`concat.nf.test` uses the text-only `concat.nf` harness to check the actual
concat script's lexical filename ordering. Shell fakes live in `bin/`;
interval preparation reuses PREPARE_INTERVALS and its GATK fake rather than
duplicating the harness.

See [the shared run instructions](../../prepare_intervals/tests/README.md)
for the Pixi-pinned nf-test command and the tests' limitations.
