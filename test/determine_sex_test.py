import unittest
import sys
sys.path.append("../")


from src.scims.determine_sex import *

tempFilesList = []

class TestDetermineSex(unittest.TestCase):

	def test_tempFilesList(self):
		with self.assertRaises(TypeError):
			tempFilesList(2)
			tempFilesList([])

	def test_alignWithBwa(self):
		"""alignWithBwa(index, forwardReads, reverseReads, threads)"""
		with self.assertRaises(TypeError):
			alignWithBwa(4, "forward", "reverse", 2)
			alignWithBwa("index", 1, "reverse", 2 )
			alignWithBwa("index", [], "reverse", 2)

	def test_calculateStats(self):
		self.assertEqual(calculateStats(2313, 1219)[0], 0.345)
		self.assertEqual(calculateStats(2313, 1219)[1], 0.016)


if __name__ == '__main__':
	test = unittest.TestDetermineSex()
	test.test_tempFilesList()