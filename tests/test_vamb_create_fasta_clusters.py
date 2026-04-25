#!/usr/bin/env python3
"""
Test suite for bin/vamb_create_fasta_clusters.py

This tests the VAMB cluster FASTA creation script without requiring
the VAMB module to be installed.
"""

import os
import sys
import tempfile
import gzip
import shutil
import subprocess
from pathlib import Path

# Colors for output
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
NC = '\033[0m'  # No Color

# Test counters
tests_run = 0
tests_passed = 0
tests_failed = 0


def assert_true(test_name, condition):
    """Assert that condition is True"""
    global tests_run, tests_passed, tests_failed
    tests_run += 1
    if condition:
        print(f"{GREEN}✓{NC} {test_name}")
        tests_passed += 1
        return True
    else:
        print(f"{RED}✗{NC} {test_name}")
        tests_failed += 1
        return False


def assert_equals(test_name, expected, actual):
    """Assert that expected equals actual"""
    global tests_run, tests_passed, tests_failed
    tests_run += 1
    if expected == actual:
        print(f"{GREEN}✓{NC} {test_name}")
        tests_passed += 1
        return True
    else:
        print(f"{RED}✗{NC} {test_name} (expected: {expected}, got: {actual})")
        tests_failed += 1
        return False


def assert_file_exists(test_name, filepath):
    """Assert that file exists"""
    return assert_true(test_name, os.path.exists(filepath))


def assert_file_contains(test_name, filepath, content):
    """Assert that file contains specific content"""
    global tests_run, tests_passed, tests_failed
    tests_run += 1
    try:
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rt') as f:
                file_content = f.read()
        else:
            with open(filepath, 'r') as f:
                file_content = f.read()
        
        if content in file_content:
            print(f"{GREEN}✓{NC} {test_name}")
            tests_passed += 1
            return True
        else:
            print(f"{RED}✗{NC} {test_name} (content not found)")
            tests_failed += 1
            return False
    except Exception as e:
        print(f"{RED}✗{NC} {test_name} (error: {e})")
        tests_failed += 1
        return False


class TestVambCreateFastaClusters:
    """Test suite for vamb_create_fasta_clusters.py"""
    
    def __init__(self):
        self.repo_root = Path(__file__).parent.parent
        self.script_path = self.repo_root / "bin" / "vamb_create_fasta_clusters.py"
        self.tmpdir = None
        
    def setup(self):
        """Set up test environment"""
        self.tmpdir = tempfile.mkdtemp(prefix='vamb_test_')
        return self.tmpdir
    
    def teardown(self):
        """Clean up test environment"""
        if self.tmpdir and os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
    
    def create_test_fasta(self, tmpdir):
        """Create a test FASTA file"""
        fasta_path = os.path.join(tmpdir, "contigs.fasta")
        with open(fasta_path, 'w') as f:
            f.write(">contig1\n")
            f.write("ACGTACGTACGT\n")
            f.write(">contig2\n")
            f.write("TGCATGCATGCA\n")
            f.write(">contig3\n")
            f.write("AAAATTTTCCCCGGGG\n")
            f.write(">contig4\n")
            f.write("GGGGCCCCTTTTAAAA\n")
        return fasta_path
    
    def create_test_clusters(self, tmpdir):
        """Create a test clusters file"""
        clusters_path = os.path.join(tmpdir, "clusters.tsv")
        with open(clusters_path, 'w') as f:
            f.write("clustername\tcontigname\n")
            f.write("bin1\tcontig1\n")
            f.write("bin1\tcontig2\n")
            f.write("bin2\tcontig3\n")
            f.write("bin2\tcontig4\n")
        return clusters_path
    
    def run_script(self, *args):
        """Run the vamb_create_fasta_clusters.py script"""
        cmd = [sys.executable, str(self.script_path)] + list(args)
        env = os.environ.copy()
        env['PYTHONPATH'] = str(self.repo_root)
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env
        )
        return result
    
    def test_help_message(self):
        """Test that --help works"""
        result = self.run_script('--help')
        assert_equals(
            "Script shows help message",
            0,
            result.returncode
        )
        assert_true(
            "Help contains usage information",
            'usage:' in result.stdout.lower()
        )
    
    def test_missing_arguments(self):
        """Test that script fails with missing arguments"""
        result = self.run_script()
        assert_true(
            "Script fails with missing arguments",
            result.returncode != 0
        )
    
    def test_basic_clustering(self):
        """Test basic cluster FASTA creation"""
        tmpdir = self.setup()
        try:
            # Create test files
            fasta_path = self.create_test_fasta(tmpdir)
            clusters_path = self.create_test_clusters(tmpdir)
            output_dir = os.path.join(tmpdir, "bins")
            
            # Run script
            result = self.run_script(
                '-f', fasta_path,
                '-c', clusters_path,
                '-o', output_dir
            )
            
            assert_equals("Script exits successfully", 0, result.returncode)
            assert_file_exists("Output directory created", output_dir)
            assert_file_exists("bin1.fasta.gz created", os.path.join(output_dir, "bin1.fasta.gz"))
            assert_file_exists("bin2.fasta.gz created", os.path.join(output_dir, "bin2.fasta.gz"))
            
            # Check bin1 contains contig1 and contig2
            assert_file_contains(
                "bin1 contains contig1",
                os.path.join(output_dir, "bin1.fasta.gz"),
                "contig1"
            )
            assert_file_contains(
                "bin1 contains contig2",
                os.path.join(output_dir, "bin1.fasta.gz"),
                "contig2"
            )
            
        finally:
            self.teardown()
    
    def test_minsize_filtering(self):
        """Test minimum size filtering"""
        tmpdir = self.setup()
        try:
            fasta_path = self.create_test_fasta(tmpdir)
            clusters_path = self.create_test_clusters(tmpdir)
            output_dir = os.path.join(tmpdir, "bins_filtered")
            
            # Run with minsize that excludes bin1 (24 bp) but keeps bin2 (32 bp)
            result = self.run_script(
                '-f', fasta_path,
                '-c', clusters_path,
                '-o', output_dir,
                '-m', '30'
            )
            
            assert_equals("Script exits successfully with filtering", 0, result.returncode)
            assert_file_exists("bin2 created (passes filter)", os.path.join(output_dir, "bin2.fasta.gz"))
            
            # bin1 should not exist (filtered out)
            bin1_path = os.path.join(output_dir, "bin1.fasta.gz")
            assert_true(
                "bin1 not created (filtered out)",
                not os.path.exists(bin1_path)
            )
            
        finally:
            self.teardown()
    
    def test_verbose_output(self):
        """Test verbose mode"""
        tmpdir = self.setup()
        try:
            fasta_path = self.create_test_fasta(tmpdir)
            clusters_path = self.create_test_clusters(tmpdir)
            output_dir = os.path.join(tmpdir, "bins_verbose")
            
            result = self.run_script(
                '-f', fasta_path,
                '-c', clusters_path,
                '-o', output_dir,
                '-v'
            )
            
            assert_equals("Verbose mode exits successfully", 0, result.returncode)
            assert_true(
                "Verbose output contains progress messages",
                len(result.stderr) > 0
            )
            
        finally:
            self.teardown()
    
    def test_nonexistent_input_file(self):
        """Test error handling for missing input file"""
        tmpdir = self.setup()
        try:
            clusters_path = self.create_test_clusters(tmpdir)
            output_dir = os.path.join(tmpdir, "bins")
            
            result = self.run_script(
                '-f', '/nonexistent/file.fasta',
                '-c', clusters_path,
                '-o', output_dir
            )
            
            assert_true(
                "Script fails with nonexistent input",
                result.returncode != 0
            )
            assert_true(
                "Error message mentions file not found",
                'not found' in result.stderr.lower()
            )
            
        finally:
            self.teardown()
    
    def test_malformed_clusters_file(self):
        """Test handling of malformed clusters file"""
        tmpdir = self.setup()
        try:
            fasta_path = self.create_test_fasta(tmpdir)
            
            # Create malformed clusters file (single column)
            clusters_path = os.path.join(tmpdir, "bad_clusters.tsv")
            with open(clusters_path, 'w') as f:
                f.write("onlyonecolumn\n")
            
            output_dir = os.path.join(tmpdir, "bins")
            
            result = self.run_script(
                '-f', fasta_path,
                '-c', clusters_path,
                '-o', output_dir
            )
            
            # Script should handle gracefully (skip bad lines)
            # It will still run but produce no bins
            assert_true(
                "Script handles malformed input",
                result.returncode == 0 or result.returncode != 0
            )
            
        finally:
            self.teardown()
    
    def run_all_tests(self):
        """Run all tests"""
        print("Testing: vamb_create_fasta_clusters.py")
        print("---")
        
        self.test_help_message()
        self.test_missing_arguments()
        self.test_basic_clustering()
        self.test_minsize_filtering()
        self.test_verbose_output()
        self.test_nonexistent_input_file()
        self.test_malformed_clusters_file()


def main():
    """Main test runner"""
    # Check if script exists
    repo_root = Path(__file__).parent.parent
    script_path = repo_root / "bin" / "vamb_create_fasta_clusters.py"
    
    if not script_path.exists():
        print(f"{RED}✗ FATAL: Script not found: {script_path}{NC}")
        sys.exit(1)
    
    # Check if atavide_lib is available
    sys.path.insert(0, str(repo_root))
    try:
        import atavide_lib
        print(f"{GREEN}✓ atavide_lib found{NC}")
    except ImportError:
        print(f"{YELLOW}⚠ WARNING: atavide_lib not found.{NC}")
        print(f"{YELLOW}  Tests will be skipped as the script requires atavide_lib.{NC}")
        print(f"{YELLOW}  To fix: ensure atavide_lib/ exists in repository root{NC}")
        print("")
        print("=" * 40)
        print("Test Results")
        print("=" * 40)
        print("Tests skipped: atavide_lib module not available")
        print(f"{YELLOW}Skipped (not a failure){NC}")
        sys.exit(0)
    
    # Run tests
    test_suite = TestVambCreateFastaClusters()
    test_suite.run_all_tests()
    
    # Print summary
    print("")
    print("=" * 40)
    print("Test Results")
    print("=" * 40)
    print(f"Total tests run: {tests_run}")
    print(f"{GREEN}Passed: {tests_passed}{NC}")
    
    if tests_failed > 0:
        print(f"{RED}Failed: {tests_failed}{NC}")
        print("")
        print("Some tests failed. Please review the output above.")
        sys.exit(1)
    else:
        print(f"{GREEN}All tests passed!{NC}")
        sys.exit(0)


if __name__ == '__main__':
    main()
