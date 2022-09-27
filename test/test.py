import glob
import os.path
import unittest

from bin.parse_chromosomes import parse_chromosomes


expected_pattern = {
    "pviv.gff3": "PvP01_([0-9][0-9])_v2",
    "padl.gff3": "PADLG01_([0-9][0-9])",
    "pfal.gff3": "Pf3D7_([0-9][0-9])_v3",
    "tgon.gff3": "TGME49_chr(.*)",
}

class TestParseChromosomes(unittest.TestCase):
    def setUp(self):
        self.test_files = glob.glob("test/*.gff3")
        self.verification_errors = []
    
    def tearDown(self):
        self.assertEqual([], self.verification_errors)
    
    def test_parse_examples(self):
        for tf in self.test_files:
            p = parse_chromosomes(tf)
            try:
                self.assertEqual(len(p.split("\t")), 2)
                self.assertEqual(p.split("\t")[1], expected_pattern[os.path.basename(tf)])
            except AssertionError as e:
                self.verification_errors.append(str(e))


if __name__ == "__main__":
    unittest.main()
