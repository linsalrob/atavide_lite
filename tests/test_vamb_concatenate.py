"""Tests for the VAMB sample mapping written by vamb_concatenate.py."""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "bin" / "vamb_concatenate.py"

# A stand-in for the vamb package (which needs ROCm PyTorch), renaming contigs the way
# vamb.vambtools.concatenate_fasta_ios does in vamb 5.0.4: S{reader_no + 1}C{identifier}
FAKE_VAMB = '''
from . import vambtools
'''
FAKE_VAMBTOOLS = '''
import gzip

def _records(path):
    opener = gzip.open if open(path, "rb").read(2) == b"\\x1f\\x8b" else open
    name, seq = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)

def concatenate_fasta(outfile, inpaths, minlength=2000, rename=True):
    for reader_no, path in enumerate(inpaths):
        for name, seq in _records(path):
            if len(seq) < minlength:
                continue
            if rename:
                name = f"S{reader_no + 1}C{name}"
            print(f">{name}\\n{seq}", file=outfile)
'''


def write_fasta(path, records, base="A"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, length in records:
            handle.write(f">{name} flag=1 len={length}\n{base * length}\n")


def read_mapping(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class VambConcatenateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        fake = self.root / "fake" / "vamb"
        fake.mkdir(parents=True)
        (fake / "__init__.py").write_text(FAKE_VAMB)
        (fake / "vambtools.py").write_text(FAKE_VAMBTOOLS)
        self.env = dict(os.environ, PYTHONPATH=str(self.root / "fake"))
        # MEGAHIT reuses contig names across samples (sample_a and sample_c both start with
        # k141_1), so only the sequence can tell a reordered input list apart. sample_b's
        # first contig is too short, so the check must skip to the next one.
        write_fasta(self.root / "megahit" / "sample_a" / "final.contigs.fa", [("k141_1", 3000), ("k141_2", 2500)])
        write_fasta(self.root / "megahit" / "sample_b" / "final.contigs.fa", [("k141_1", 100), ("k141_7", 2200)], base="C")
        write_fasta(self.root / "megahit" / "sample_c" / "output" / "final.contigs.fa", [("k141_1", 4000)], base="G")
        self.inputs = [
            str(self.root / "megahit" / "sample_a" / "final.contigs.fa"),
            str(self.root / "megahit" / "sample_b" / "final.contigs.fa"),
            str(self.root / "megahit" / "sample_c" / "output" / "final.contigs.fa"),
        ]
        (self.root / "vamb").mkdir()
        self.contigs = self.root / "vamb" / "contigs.fna.gz"
        self.mapping = self.root / "vamb" / "sample_mappings.tsv"

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, *arguments):
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            env=self.env, capture_output=True, text=True,
        )

    def test_concatenate_writes_mapping_matching_contig_prefixes(self):
        result = self.run_script("-m", "2000", str(self.contigs), *self.inputs)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = read_mapping(self.mapping)
        self.assertEqual([(r["vamb_sample"], r["sample"]) for r in rows],
                         [("S1", "sample_a"), ("S2", "sample_b"), ("S3", "sample_c")])
        self.assertEqual([r["path"] for r in rows], self.inputs)
        with gzip.open(self.contigs, "rt") as handle:
            headers = [line[1:].strip() for line in handle if line.startswith(">")]
        self.assertEqual(headers, ["S1Ck141_1", "S1Ck141_2", "S2Ck141_7", "S3Ck141_1"])

    def test_mapping_only_recreates_mapping_for_existing_output(self):
        self.assertEqual(self.run_script(str(self.contigs), *self.inputs).returncode, 0)
        expected = self.mapping.read_text()
        self.mapping.unlink()
        result = self.run_script("--mapping-only", str(self.contigs), *self.inputs)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.mapping.read_text(), expected)

    def test_mapping_only_rejects_inputs_in_a_different_order(self):
        self.assertEqual(self.run_script(str(self.contigs), *self.inputs).returncode, 0)
        self.mapping.unlink()
        # Swapping sample_a and sample_c keeps every expected name (S1Ck141_1, S3Ck141_1)
        # present, so this is caught only by comparing sequences
        swapped = [self.inputs[2], self.inputs[1], self.inputs[0]]
        result = self.run_script("--mapping-only", str(self.contigs), *swapped)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not in the order", result.stderr)
        self.assertFalse(self.mapping.exists())

    def test_mapping_only_does_not_need_vamb(self):
        self.assertEqual(self.run_script(str(self.contigs), *self.inputs).returncode, 0)
        self.mapping.unlink()
        self.env["PYTHONPATH"] = ""
        result = self.run_script("--mapping-only", str(self.contigs), *self.inputs)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_refuses_to_overwrite_existing_mapping(self):
        self.mapping.write_text("old\n")
        result = self.run_script(str(self.contigs), *self.inputs)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.contigs.exists())
        self.assertEqual(self.mapping.read_text(), "old\n")

    def test_keepnames_writes_no_mapping(self):
        result = self.run_script("--keepnames", str(self.contigs), *self.inputs)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.mapping.exists())

    def test_sample_name_from_file_name(self):
        write_fasta(self.root / "assemblies" / "sample_d.fa.gz", [("c1", 3000)])
        out = self.root / "named" / "contigs.fna.gz"
        out.parent.mkdir()
        # write_fasta writes plain text; the gz suffix only exercises the name parsing
        result = self.run_script(str(out), str(self.root / "assemblies" / "sample_d.fa.gz"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(read_mapping(out.parent / "sample_mappings.tsv")[0]["sample"], "sample_d")


if __name__ == "__main__":
    unittest.main()
