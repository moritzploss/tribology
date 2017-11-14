import unittest
from tribology import dyn2kin


class TestUM(unittest.TestCase):

    def test_dyn2kin(self):
        self.assertEqual(dyn2kin(1, 0.5), 0.5)
