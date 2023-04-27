import glob
import os.path
import unittest

from bin.parse_chromosomes import ChromosomeParser


expected_pattern = {
    "pviv.gff3": "PvP01_([0-9][0-9]){0,1}(MIT|API){0,1}_v2",
    "padl.gff3": "PADLG01_([0-9][0-9])",
    "pfal.gff3": "Pf3D7_([0-9][0-9]){0,1}(MIT|API){0,1}_v3",
    "tgon.gff3": "TGME49_([MDCLXVI]+[a-z]?)",
    "ltro.gff3": "LtrL590_([0-9][0-9])",
    "lmxm.gff3": "LmxM.([0-9][0-9])",
    "lael.gff3": "LaeL147_([0-9][0-9])",
    "larl.gff3": "LarLEM1108_([0-9][0-9])",
    "prel.gff3": "PRELSG_([0-9][0-9]){0,1}(MIT|API){0,1}_v1"
}


class TestParseChromosomes(unittest.TestCase):
    def setUp(self):
        self.test_files = glob.glob("test/*.gff3")
        self.verification_errors = []

    def tearDown(self):
        self.assertEqual([], self.verification_errors)

    def test_parse_examples(self):
        for tf in self.test_files:
            p = ChromosomeParser(tf)
            p.generate_regex()
            p.correct_prefix()
            try:
                self.assertEqual(len(str(p).split("\t")), 2)
                self.assertEqual(str(p).split("\t")[1],
                                 expected_pattern[os.path.basename(tf)])
            except AssertionError as e:
                self.verification_errors.append(str(e))


if __name__ == "__main__":
    unittest.main()
