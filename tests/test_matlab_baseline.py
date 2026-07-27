import unittest
from pathlib import Path


class MatlabBaselineArtifactsTest(unittest.TestCase):
    def test_baslab_entrypoint_exists(self):
        repo_root = Path(__file__).resolve().parents[1]
        script = repo_root / 'matlab' / 'createMinimalBaseline.m'
        self.assertTrue(script.exists(), f'Missing MATLAB baseline entrypoint: {script}')


if __name__ == '__main__':
    unittest.main()
