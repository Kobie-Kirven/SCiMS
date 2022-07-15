
import unittest
import sys

sys.path.append("../")


from src.scims.determine_sex import *


class TestDetermineSex(unittest.TestCase):

    def test_get_align_handle(self):
        type_file = "<class 'pysam.libcalignmentfile.AlignmentFile'>"
        self.assertEqual(str(get_align_handle("../test_data/male_test.sam").__class__), type_file)
        self.assertEqual(str(get_align_handle("../test_data/male_test_1.bam").__class__), type_file)
        with self.assertRaises(WrongFile):
            get_align_handle("../test_data/scaffolds.txt")
    
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

if __name__ == "__main__":
    unittest.main()
