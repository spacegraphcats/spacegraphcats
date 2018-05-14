import unittest

from .frontier_search import find_shadow


class SearchTest(unittest.TestCase):
    def test_find_shadow(self):
        dag = {0: [1, 2], 1: [3, 4, 6], 2: [5, 6], 3: [], 4: [], 5: [], 6: []}

        self.assertEqual(find_shadow([0], dag), set([3, 4, 5, 6]))

if __name__ == '__main__':
    unittest.main()
