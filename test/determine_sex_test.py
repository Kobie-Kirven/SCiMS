import unittest
import sys
sys.path.append("../")
import Bio
from src.scims.determine_sex import *

tempFilesList = []

class TestDetermineSex(unittest.TestCase):

	def test_tempFilesList(self):

		with self.assertRaises(TypeError):
			tempFilesList(2)
			tempFilesList([])


if __name__ == '__main__':
	test = unittest.TestDetermineSex()
	test.test_tempFilesList()