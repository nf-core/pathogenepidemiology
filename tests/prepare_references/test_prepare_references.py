"""Run with: python -m unittest discover -s tests/prepare_references -v.

Requires Pixi's pinned Nextflow launcher and local wget, gzip, samtools and gatk.
All fixtures are synthetic; HTTP requests stay on localhost. Logs and task
artifacts are retained under /tmp/opencode for debugging.
"""

import functools
import gzip
import http.server
import json
import re
import shutil
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FASTA = ">synthetic\nACGTACGT\n"


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass


class PrepareReferencesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        missing = [
            tool
            for tool in ("pixi", "wget", "gzip", "samtools", "gatk")
            if not shutil.which(tool)
        ]
        if missing:
            raise unittest.SkipTest(
                f"Missing local test dependencies: {', '.join(missing)}"
            )
        cls.artifacts = Path(
            tempfile.mkdtemp(prefix="prepare-references-", dir="/tmp/opencode")
        )
        cls.fixtures = cls.artifacts / "fixtures"
        cls.fixtures.mkdir()
        for extension in ("fa", "fasta", "fna"):
            (cls.fixtures / f"reference.{extension}").write_text(FASTA)
        with gzip.open(cls.fixtures / "reference.fna.gz", "wt") as handle:
            handle.write(FASTA)
        with gzip.open(cls.fixtures / "remote.fa.gz", "wt") as handle:
            handle.write(FASTA)
        (cls.fixtures / "alternate").mkdir()
        with gzip.open(cls.fixtures / "alternate" / "remote.fa.gz", "wt") as handle:
            handle.write(">synthetic\nACGT\n")
        (cls.fixtures / "host.fa").write_text(">host\nTTTT\n")
        with gzip.open(cls.fixtures / "host.fa.gz", "wt") as handle:
            handle.write(">host\nTTTT\n")
        cls.server = http.server.ThreadingHTTPServer(
            ("127.0.0.1", 0),
            functools.partial(QuietHandler, directory=str(cls.fixtures)),
        )
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base_url = f"http://127.0.0.1:{cls.server.server_port}"
        print(f"\nReference test artifacts: {cls.artifacts}", flush=True)

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()

    def run_reference(
        self,
        query=None,
        remote="unused.fa.gz",
        extra_args=(),
        host=None,
        host_remote="host.fa",
        run_dir=None,
    ):
        if run_dir is None:
            run_dir = Path(tempfile.mkdtemp(prefix="run-", dir=self.artifacts))
        command = [
            "pixi",
            "run",
            "nextflow",
            "-log",
            str(run_dir / "nextflow.log"),
            "-C",
            str(HERE / "nextflow.config"),
            "run",
            str(HERE / "main.nf"),
            "-ansi-log",
            "false",
            "-work-dir",
            str(run_dir / "work"),
            "-with-trace",
            str(run_dir / "trace.tsv"),
            "--outdir",
            str(run_dir / "results"),
            "--queryurl",
            f"{self.base_url}/{remote}",
            "--validate_params",
            "false",
        ]
        if host_remote is not None:
            command.extend(["--hosturl", f"{self.base_url}/{host_remote}"])
        else:
            params_file = run_dir / "no-host.json"
            params_file.write_text(json.dumps({"hosturl": None}))
            command.extend(["-params-file", str(params_file)])
        if query is not None:
            command.extend(["--query_fasta", str(query)])
        if host is not None:
            command.extend(["--host_fasta", str(host)])
        command.extend(extra_args)
        result = subprocess.run(
            command, cwd=ROOT, capture_output=True, text=True, timeout=120, check=False
        )
        output = result.stdout + result.stderr
        (run_dir / "output.log").write_text(output)
        self.last_run = run_dir
        trace_path = run_dir / "trace.tsv"
        trace = trace_path.read_text() if trace_path.exists() else ""
        return result.returncode, output, trace

    def assert_reference(
        self,
        result,
        name,
        compressed=False,
        downloaded=False,
        host_compressed=False,
        host_downloaded=True,
    ):
        code, output, trace = result
        self.assertEqual(code, 0, output)
        meta = Path(name).stem
        self.assertIn(f"QUERY\t{meta}\t{name}\t{FASTA.strip()}", output)
        self.assertIn(f"FAI\t{meta}\t{name}.fai\tsynthetic\t8\t11\t8\t9", output)
        self.assertIn("HOST\thost.fa\t>host\nTTTT", output)
        self.assertIn("SAMTOOLS_FAIDX", trace)
        self.assertIn("GATK4_CREATESEQUENCEDICTIONARY", trace)
        self.assertIn(f"DICT\t{meta}\t{meta}.dict", output)
        self.assertEqual("DOWNLOAD_HOSTFASTA" in trace, host_downloaded)
        self.assertEqual("GUNZIP_HOSTFASTA" in trace, host_compressed)
        self.assertEqual("GUNZIP_QUERYFASTA" in trace, compressed)
        self.assertEqual("DOWNLOAD_QUERYFASTA" in trace, downloaded)
        published = self.last_run / "results" / "download"
        self.assertEqual(
            (published / f"{name}.fai").read_text(), "synthetic\t8\t11\t8\t9\n"
        )
        if downloaded or compressed:
            self.assertEqual((published / name).read_text(), FASTA)
        if host_downloaded or host_compressed:
            self.assertEqual((published / "host.fa").read_text(), ">host\nTTTT\n")

    def test_both_local_plain_or_gzip_skip_both_downloads(self):
        for query_gz in (False, True):
            for host_gz in (False, True):
                with self.subTest(query_gz=query_gz, host_gz=host_gz):
                    query = "reference.fna.gz" if query_gz else "reference.fa"
                    host = "host.fa.gz" if host_gz else "host.fa"
                    self.assert_reference(
                        self.run_reference(
                            self.fixtures / query,
                            host=self.fixtures / host,
                            host_remote="does-not-exist.fa.gz",
                        ),
                        "reference.fna" if query_gz else "reference.fa",
                        compressed=query_gz,
                        host_compressed=host_gz,
                        host_downloaded=False,
                    )

    def test_download_query_with_local_host(self):
        self.assert_reference(
            self.run_reference(remote="remote.fa.gz", host=self.fixtures / "host.fa"),
            "remote.fa",
            compressed=True,
            downloaded=True,
            host_downloaded=False,
        )

    def test_absent_host_preserves_empty_value(self):
        code, output, trace = self.run_reference(
            self.fixtures / "reference.fa", host_remote=None
        )
        self.assertEqual(code, 0, output)
        self.assertIn("HOST\tabsent", output)
        self.assertNotIn("DOWNLOAD_HOSTFASTA", trace)
        self.assertNotIn("GUNZIP_HOSTFASTA", trace)

    def test_local_host_overrides_absent_url(self):
        self.assert_reference(
            self.run_reference(
                self.fixtures / "reference.fa",
                host=self.fixtures / "host.fa",
                host_remote=None,
            ),
            "reference.fa",
            host_downloaded=False,
        )

    def test_both_downloaded_gzip(self):
        self.assert_reference(
            self.run_reference(remote="remote.fa.gz", host_remote="host.fa.gz"),
            "remote.fa",
            compressed=True,
            downloaded=True,
            host_compressed=True,
        )

    def test_changed_url_with_same_basename_is_not_reused_on_resume(self):
        self.assert_reference(
            self.run_reference(remote="remote.fa.gz", host=self.fixtures / "host.fa"),
            "remote.fa",
            compressed=True,
            downloaded=True,
            host_downloaded=False,
        )
        run_dir = self.last_run
        session = re.search(
            r"Session UUID: ([\w-]+)", (run_dir / "nextflow.log").read_text()
        )
        self.assertIsNotNone(session)
        code, output, trace = self.run_reference(
            remote="alternate/remote.fa.gz",
            host=self.fixtures / "host.fa",
            run_dir=run_dir,
            extra_args=("-resume", session.group(1)),
        )
        self.assertEqual(code, 0, output)
        self.assertIn("QUERY\tremote\tremote.fa\t>synthetic\nACGT", output)
        self.assertIn("FAI\tremote\tremote.fa.fai\tsynthetic\t4\t11\t4\t5", output)
        self.assertIn("DOWNLOAD_QUERYFASTA", trace)
        self.assertNotIn("CACHED", trace)
        self.assertEqual(
            (run_dir / "results" / "download" / "remote.fa").read_text(),
            ">synthetic\nACGT\n",
        )

    def test_plain_downloads_bypass_gunzip(self):
        self.assert_reference(
            self.run_reference(remote="reference.fa"),
            "reference.fa",
            downloaded=True,
        )

    def test_host_invalid_inputs_fail_without_validation(self):
        directory = self.fixtures / "host-directory.fa"
        directory.mkdir()
        for host in (
            self.fixtures / "missing-host.fa",
            directory,
            self.fixtures / "*.fa",
            self.fixtures / "host.txt",
        ):
            with self.subTest(host=host):
                code, output, _ = self.run_reference(
                    self.fixtures / "reference.fa", host=host
                )
                self.assertNotEqual(code, 0, output)
                self.assertRegex(
                    output, r"host_fasta|No such file or directory|does not exist"
                )

    def test_host_list_fails_without_validation(self):
        params_file = self.artifacts / "multiple-hosts.json"
        params_file.write_text(
            json.dumps({"host_fasta": [str(self.fixtures / "host.fa")]})
        )
        code, output, _ = self.run_reference(
            self.fixtures / "reference.fa",
            extra_args=("-params-file", str(params_file)),
        )
        self.assertNotEqual(code, 0, output)
        self.assertIn("--host_fasta must be a single local FASTA file", output)

    def test_plain_local_formats_override_url_and_generate_index(self):
        for extension in ("fa", "fasta", "fna"):
            with self.subTest(extension=extension):
                name = f"reference.{extension}"
                self.assert_reference(self.run_reference(self.fixtures / name), name)

    def test_gzip_local_overrides_url_and_generates_index(self):
        self.assert_reference(
            self.run_reference(self.fixtures / "reference.fna.gz"),
            "reference.fna",
            compressed=True,
        )

    def test_url_fallback_still_downloads_and_decompresses(self):
        self.assert_reference(
            self.run_reference(remote="remote.fa.gz"),
            "remote.fa",
            compressed=True,
            downloaded=True,
        )

    def test_missing_path_fails_without_schema_validation(self):
        code, output, _trace = self.run_reference(self.fixtures / "missing.fa")
        self.assertNotEqual(code, 0, output)
        self.assertIn("missing.fa", output)
        self.assertRegex(output, r"No such file or directory|does not exist")

    def test_directory_fails_without_schema_validation(self):
        directory = self.fixtures / "directory.fa"
        directory.mkdir()
        code, output, _trace = self.run_reference(directory)
        self.assertNotEqual(code, 0, output)
        self.assertIn("--query_fasta must be a single file", output)

    def test_wildcard_does_not_expand_to_multiple_files(self):
        code, output, _trace = self.run_reference(self.fixtures / "*.fa")
        self.assertNotEqual(code, 0, output)
        self.assertRegex(output, r"No such file or directory|does not exist")

    def test_invalid_extension_fails_without_schema_validation(self):
        invalid = self.fixtures / "reference.txt"
        invalid.write_text(FASTA)
        code, output, _trace = self.run_reference(invalid)
        self.assertNotEqual(code, 0, output)
        self.assertIn("--query_fasta must be a single local FASTA file", output)

    def test_list_fails_without_schema_validation(self):
        params_file = self.artifacts / "multiple-references.json"
        params_file.write_text(
            json.dumps(
                {
                    "query_fasta": [
                        str(self.fixtures / "reference.fa"),
                        str(self.fixtures / "reference.fna"),
                    ]
                }
            )
        )
        code, output, _trace = self.run_reference(
            extra_args=("-params-file", str(params_file))
        )
        self.assertNotEqual(code, 0, output)
        self.assertIn("--query_fasta must be a single local FASTA file", output)

    def test_unsafe_local_basenames_fail_before_tasks(self):
        for option in ("query_fasta", "host_fasta"):
            for stem in ("unsafe name", "unsafe;name", "unsafe$(name)", "-unsafe"):
                for extension in ("fa", "fa.gz"):
                    with self.subTest(option=option, stem=stem, extension=extension):
                        unsafe = self.fixtures / f"{stem}.{extension}"
                        if extension.endswith(".gz"):
                            with gzip.open(unsafe, "wt") as handle:
                                handle.write(FASTA)
                        else:
                            unsafe.write_text(FASTA)
                        query = (
                            unsafe
                            if option == "query_fasta"
                            else self.fixtures / "reference.fa"
                        )
                        host = (
                            unsafe
                            if option == "host_fasta"
                            else self.fixtures / "host.fa"
                        )
                        code, output, trace = self.run_reference(query, host=host)
                        self.assertNotEqual(code, 0, output)
                        self.assertIn(f"--{option} basename must start", output)
                        self.assertIn(unsafe.name, output)
                        self.assertNotIn("DOWNLOAD_", trace)
                        self.assertNotIn("GUNZIP_", trace)
                        self.assertNotIn("SAMTOOLS_FAIDX", trace)
                        self.assertNotIn("GATK4_CREATESEQUENCEDICTIONARY", trace)

    def test_schema_declares_optional_existing_fasta_file(self):
        schema = json.loads((ROOT / "nextflow_schema.json").read_text())
        options = schema["$defs"]["pipeline_input_options"]
        query = options["properties"]["query_fasta"]
        host = options["properties"]["host_fasta"]
        self.assertNotIn("host_fasta", options.get("required", []))
        for field in ("type", "format", "exists", "pattern"):
            self.assertEqual(query[field], host[field])
        self.assertNotIn("query_fasta", options.get("required", []))
        self.assertEqual(query["type"], "string")
        self.assertEqual(query["format"], "file-path")
        self.assertTrue(query["exists"])
        for extension in ("fa", "fasta", "fna", "fa.gz", "fasta.gz", "fna.gz"):
            self.assertIsNotNone(
                re.fullmatch(query["pattern"], f"reference.{extension}")
            )
        self.assertIsNone(re.fullmatch(query["pattern"], "reference.txt"))


if __name__ == "__main__":
    unittest.main()
