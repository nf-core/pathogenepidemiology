# nf-core/pathogenepidemiology: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

In the GATK subworkflow, chromosome VCFs larger than 8000 bytes pass through VQSR.
VCFs of 8000 bytes or smaller bypass VQSR and contribute **unrecalibrated calls**
to the cohort alongside completed VQSR outputs; small chromosomes are not dropped.
VCFs remain paired with their indexes. Each contig has one calling BED and one
cohort VCF: `chr1.vcf.gz` (or `chr1.snp.recal.vcf.gz` after VQSR).
Per-sample GVCFs are named `sample.chr1.g.vcf.gz`; the final cohort name is unchanged.
Duplicate contigs across interval files cause a workflow error before conversion.
`GATK_MOI` collects the indexed VQSR and bypass outputs directly into one cohort.
The nf-core `BCFTOOLS_CONCAT` wrapper sorts input VCF filenames lexicographically
(`chr1`, `chr10`, `chr2`) before invoking `bcftools concat`; the executable does
not perform this sorting. The pipeline assumes this filename order is the intended
concatenation order, not reference-index or natural chromosome order. Original VCF
and index filenames are retained without numbered staging aliases.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC for non-ONT samples
- [NanoPlot and Filtlong](#nanoplot-and-filtlong) - ONT raw read QC and length filtering
- [BBDuk](#bbduk) - Non-ONT adapter trimming and contaminant filtering
- [ONT alignment](#ont-alignment) - Competitive mapping and selection of primary parasite alignments
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### BBDuk

For non-ONT samples, the local `bbduk_custom` module runs as two distinct tasks,
`BBDUK_TRIM` followed by `BBDUK_FILTER`.

- `trimmed/`: adapter-trimmed `*_trimmed.fastq.gz`, `*.trimmed.bbduk.log`,
  `*.trimmed.stats.txt` and `*.trimmed.discarded.fastq.gz`. These rejects failed
  trimming/length checks; they are **not** contaminant matches.
- `filtered/`: final `*_filtered.fastq.gz`, `*.filtered.bbduk.log`,
  `*.filtered.stats.txt` and `*.filtered.discarded.fastq.gz`. Filtering removes whole
  reads matching PhiX or the supplied host/contaminant reference, rather than
  trimming them. For paired-end data the default BBDuk pair-removal policy is
  retained: a matching mate rejects the pair, so a rejected mate need not itself match.

The preprocessing `bbduk_stats` and `bbduk_logs` channels contain reports from
all BBDuk stages with distinct filenames. `bbduk_dropped` contains contaminant
filtering rejects only. Trimming rejects are still published separately in `trimmed/`.

### NanoPlot and Filtlong

- `nanoplot/<run_accession>/`: Sample-prefixed NanoPlot reports, plots and
  `*NanoStats.txt` from raw ONT reads. MultiQC includes these raw-read statistics.
  MultiQC uses the original run ID to join NanoStats and ONT SAMtools metrics,
  subject to MultiQC's filename cleaning. IDs containing recognised filename
  tokens and custom sample-renaming or filename-cleaning settings can affect
  the displayed names and joins.
- `filtlong/`: `<run_accession>.filtlong.fastq.gz` and
  `<run_accession>.filtlong.log`. The length-filtered FASTQ is the input to
  competitive alignment and has not been filtered against a host reference.
  MultiQC does not summarise the default length-only Filtlong run: its parser
  requires a target base count in the log. The log remains available in `filtlong/`.

### ONT alignment

- `aligned_l/<run_accession>.bam` and `.bam.bai`: coordinate-sorted, indexed BAM
  of mapped primary parasite alignments overlapping the core BED. Overlapping
  alignments are retained untrimmed; off-core alignments are excluded. The
  BAM header uses the parasite reference when no host is supplied. With a host,
  it retains original contig IDs from the combined reference, including
  host IDs, but no alignments to host contigs are retained. An empty result
  still has a BAM and index.

MultiQC's ONT SAMtools alignment statistics are computed from the retained BAM,
not the combined alignment or the raw FASTQ. A selected alignment does not rule
out host sequence in the read. No host-depleted FASTQ is produced. The competitive
BAM is retained in the task work directory and is not published. The combined
reference and minimap2 index are generated in task work directories and may also
be published by process configuration.
`FIND_CONCATENATE` orders the original
parasite and host FASTAs by staged basename, so the combined reference may place
host records before parasite records.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
