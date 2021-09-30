import unittest
import sys
import pandas as pd

sys.path.append("../")


from src.scims.determine_sex import *

tempFilesList = []


class TestDetermineSex(unittest.TestCase):
    def test_check_index_files(self):
        self.assertEqual(
            check_index_files(
                ["amb", "ann", "bwt", "pac", "sa"],
                "/Users/kobiekirven/Desktop/scimsTest/trim",
            ),
            True,
        )
        self.assertEqual(
            check_index_files(
                ["amb", "ann", "bwt", "pac", "sa"],
                "/Users/kobiekirven/Desktop/scimsTest/trim1",
            ),
            False,
        )
        os.chdir("/Users/kobiekirven/Desktop/scimsTest/")
        self.assertEqual(
            check_index_files(["amb", "ann", "bwt", "pac", "sa"], "trim"), True
        )

    def test_not_close_to_max(self):
        my_list = [{"read_name": "hello", "alignment_score": 1}, {"read_name": "hell", "alignment_score": 4},
                   {"read_name": "hello", "alignment_score": 600},
                   {"read_name": "hello", "alignment_score": 7}, {"read_name": "hell", "alignment_score": 5},
                   {"read_name": "hello", "alignment_score": 8}]
        df = pd.DataFrame(my_list)
        print(not_close_to_max(df, 0.025))

    def test_check_bwa_index(self):
        with self.assertRaises(Exception):
            check_bwa_index("trim")
        self.assertEqual(check_bwa_index("/Users/kobiekirven/Desktop/scimsTest/trim"), True)

    def test_temp_files_list(self):
        with self.assertRaises(TypeError):
            temp_files_list(2)
            temp_files_list([])

    def test_align_with_bwa(self):
        """align_with_bwa(index, forward_reads, reverse_reads, threads)"""
        with self.assertRaises(TypeError):
            align_with_bwa(4, "forward", "reverse", 2)
            align_with_bwa("index", 1, "reverse", 2)
            align_with_bwa("index", [], "reverse", 2)

    def test_calculate_stats(self):
        self.assertEqual(calculate_stats(2313, 1219)[0], 0.345)
        self.assertEqual(calculate_stats(2313, 1219)[1], 0.016)


if __name__ == "__main__":
    test = unittest.TestDetermineSex()
    test.test_check_index_files
    test.test_not_close_to_max()
    test.test_tempFilesList()
    test.test_check_bwa_inde()
