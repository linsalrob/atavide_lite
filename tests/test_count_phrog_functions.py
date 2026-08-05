"""Focused tests for read-level PHROG marker counting."""

import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "bin" / "count_phrog_functions.py"

ANNOTATIONS = [
    ("1", "integrase", "integration and excision"),
    ("phrog_2", "Cox-like excisionase and repressor", "transcription regulation"),
    ("3", "helicase subunit of the DNA excision repair complex", "DNA metabolism"),
    ("4", "transcriptional regulator", "transcription regulation"),
    ("5", "major capsid protein", "head and packaging"),
    ("6", "spore coat protein", "other"),
    ("7", "minor coat protein", "head and packaging"),
    ("8", "portal protein", "connector"),
    ("9", "terminase large subunit", "head and packaging"),
    ("10", "tail associated lysin", "tail"),
    ("11", "holin/anti-holin", "lysis"),
    ("12", "endolysin", "lysis"),
    ("13", "hemolysin", "other"),
    ("14", "", "unknown function"),
    ("15", "recombination directionality factor", "integration and excision"),
]


def report_row(count, phrog="", annotation="", category="", color=""):
    return ["target", str(count), "x", "x", "x", "x", "x", "x", phrog, color, annotation, category]


def read_matrix(path):
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    samples = rows[0][1:]
    return {row[0]: dict(zip(samples, map(float, row[1:]))) for row in rows[1:]}


class CountPhrogFunctionsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.database = self.root / "phrog.tsv.gz"
        with gzip.open(self.database, "wt", encoding="utf-8") as output:
            output.write("phrog\tcolor\tannot\tcategory\n")
            for phrog, annotation, category in ANNOTATIONS:
                output.write(f"{phrog}\tblue\t{annotation}\t{category}\n")

        rows = [
            report_row(2, "phrog_1", "integrase", "integration and excision", "blue"),
            report_row(3, "phrog_2", "Cox-like excisionase and repressor", "transcription regulation", "blue"),
            report_row(5, "phrog_3", ANNOTATIONS[2][1], "DNA metabolism", "blue"),
            report_row(7, "phrog_4", "transcriptional regulator", "transcription regulation", "blue"),
            report_row(11, "phrog_5", "major capsid protein", "head and packaging", "blue"),
            report_row(13, "phrog_6", "spore coat protein", "other", "blue"),
            report_row(17, "phrog_7", "minor coat protein", "head and packaging", "blue"),
            report_row(19, "phrog_8", "portal protein", "connector", "blue"),
            report_row(23, "phrog_9", "terminase large subunit", "head and packaging", "blue"),
            report_row(29, "phrog_10", "tail associated lysin", "tail", "blue"),
            report_row(31, "phrog_11", "holin/anti-holin", "lysis", "blue"),
            report_row(37, "phrog_12", "endolysin", "lysis", "blue"),
            report_row(41, "phrog_13", "hemolysin", "other", "blue"),
            report_row(43, "phrog_14", "", "unknown function", "blue"),
            report_row(47, "phrog_15", "recombination directionality factor", "integration and excision", "blue"),
            report_row(67),
        ]
        self.write_report("sample_a", rows)
        self.write_report("sample_zero", [report_row(100)])
        self.output = self.root / "output"

    def tearDown(self):
        self.temporary.cleanup()

    def write_report(self, sample, rows):
        directory = self.root / "mmseqs" / sample
        directory.mkdir(parents=True)
        path = directory / f"{sample}_tophit_report_phrog_function.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8") as output:
            for row in rows:
                output.write("\t".join(row) + "\n")

    def run_script(self, database=None, expect_success=True):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "-d", str(self.root / "mmseqs"), "-o", str(self.output), "-n", "test", "-a", str(database or self.database)],
            text=True,
            capture_output=True,
        )
        if expect_success and result.returncode:
            self.fail(result.stderr)
        return result

    def test_outputs_counts_normalisation_exclusions_and_zero_sample(self):
        self.run_script()
        expected_legacy = {
            "test_phrog_raw.tsv", "test_phrog_norm_all.tsv", "test_phrog_norm_phrog.tsv",
            "test_annotation_raw.tsv", "test_annotation_norm_all.tsv", "test_annotation_norm_phrog.tsv",
            "test_category_raw.tsv", "test_category_norm_all.tsv", "test_category_norm_phrog.tsv",
            "test_phrog_metadata.tsv",
        }
        self.assertTrue(expected_legacy.issubset({path.name for path in self.output.iterdir()}))
        with (self.output / "test_phrog_raw.tsv").open() as handle:
            self.assertEqual(handle.readline().rstrip("\n"), "\tsample_a\tsample_zero")
            self.assertEqual(handle.readline().rstrip("\n"), "phrog_1\t2\t0")

        raw = read_matrix(self.output / "test_marker_raw.tsv")
        expected = {
            "all_phrog_hits": 328, "known_function_phrog_hits": 285,
            "unknown_function_phrog_hits": 43, "integrase": 2, "excisionase": 50,
            "integrase_plus_excisionase": 52, "lysogeny_associated_regulators": 3,
            "major_capsid": 11, "portal": 19, "terminase": 23, "tail": 29,
            "holin": 31, "endolysin": 66, "holin_plus_endolysin": 97,
        }
        for marker, count in expected.items():
            self.assertEqual(raw[marker]["sample_a"], count)
            self.assertEqual(raw[marker]["sample_zero"], 0)

        norm_all = read_matrix(self.output / "test_marker_norm_all.tsv")
        norm_phrog = read_matrix(self.output / "test_marker_norm_phrog.tsv")
        self.assertAlmostEqual(
            norm_all["all_phrog_hits"]["sample_a"], 328 * 1_000_000 / 395, places=5
        )
        self.assertEqual(norm_phrog["all_phrog_hits"]["sample_a"], 1_000_000)
        self.assertEqual(norm_phrog["all_phrog_hits"]["sample_zero"], 0)

        with (self.output / "test_phrog_count_summary.tsv").open() as handle:
            summary = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(float(summary["sample_a"]["known_function_phrog_hits"]) + float(summary["sample_a"]["unknown_function_phrog_hits"]), float(summary["sample_a"]["all_phrog_hits"]))

        mapping = (self.output / "test_phrog_marker_mapping.tsv").read_text()
        self.assertIn("phrog_2\tCox-like excisionase and repressor\ttranscription regulation\texcisionase", mapping)
        self.assertIn("phrog_2\tCox-like excisionase and repressor\ttranscription regulation\tlysogeny_associated_regulators", mapping)
        excluded = (
            "phrog_3\t" + ANNOTATIONS[2][1] + "\tDNA metabolism\texcisionase",
            "phrog_4\ttranscriptional regulator\ttranscription regulation\tlysogeny_associated_regulators",
            "phrog_6\tspore coat protein\tother\tmajor_capsid",
            "phrog_7\tminor coat protein\thead and packaging\tmajor_capsid",
            "phrog_13\themolysin\tother\tholin",
            "phrog_13\themolysin\tother\tendolysin",
        )
        for relationship in excluded:
            self.assertNotIn(relationship, mapping)

    def test_uncompressed_annotation_input(self):
        plain = self.root / "phrog.tsv"
        with gzip.open(self.database, "rt", encoding="utf-8") as source:
            plain.write_text(source.read(), encoding="utf-8")
        self.run_script(plain)

    def test_malformed_and_conflicting_annotation_records_fail_clearly(self):
        malformed = self.root / "malformed.tsv"
        malformed.write_text("phrog\tannot\n1\tintegrase\n", encoding="utf-8")
        result = self.run_script(malformed, expect_success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("must contain tab-separated columns", result.stderr)

        conflicting = self.root / "conflicting.tsv"
        conflicting.write_text("phrog\tannot\tcategory\n1\tintegrase\tintegration\nphrog_1\tportal protein\tconnector\n", encoding="utf-8")
        result = self.run_script(conflicting, expect_success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("conflicting record for PHROG 1", result.stderr)


if __name__ == "__main__":
    unittest.main()
