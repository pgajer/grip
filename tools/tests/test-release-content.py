"""Exercise repository checks in isolated source fixtures, not the working tree."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

CHECK = Path(__file__).resolve().parents[1] / 'check-release-content.R'

class ReleaseContentTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'man').mkdir()
        (self.root / 'DESCRIPTION').write_text('Package: fixture\nVersion: 0.0.1\n')
        (self.root / 'NAMESPACE').write_text('export(example)\n')
        (self.root / 'man/example.Rd').write_text(
            '\\name{example}\n\\alias{example}\n\\title{Example}\n\\description{Example.}\n')
        subprocess.run(['git', 'init', '-q'], cwd=self.root, check=True)

    def check(self, expected=None):
        subprocess.run(['git', 'add', '.'], cwd=self.root, check=True)
        result = subprocess.run(['Rscript', str(CHECK)], cwd=self.root,
                                text=True, capture_output=True)
        if expected is None:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(expected, result.stdout + result.stderr)

    def test_documented_export_without_private_catalog(self):
        self.check()

    def test_missing_documentation_directory(self):
        shutil.rmtree(self.root / 'man')
        self.check('must both be present')

    def test_missing_export(self):
        (self.root / 'NAMESPACE').write_text('')
        self.check('must both be present')

    def test_undocumented_export(self):
        (self.root / 'NAMESPACE').write_text('export(example)\nexport(undocumented)\n')
        self.check('Public exports lack package Rd aliases: undocumented')

    def test_private_reference(self):
        (self.root / 'README.md').write_text('Private link: ' + '.codex/' + 'private/fixture\n')
        self.check('Public content refers to private working material')

    def test_private_path_in_serialized_attribute(self):
        subprocess.run(['Rscript', '-e',
            'saveRDS(structure(list(value=1), provenance="/Users/fixture/input"), "data.rds")'],
            cwd=self.root, check=True, capture_output=True)
        self.check('Personal filesystem path in serialized data')

if __name__ == '__main__':
    unittest.main()
