import utils
import unittest
import itertools

class UtilsCase(unittest.TestCase):

    def test_row_to_block(self):
        map = {10:1,
               15:2,
               30:3}

        for row in map:
            self.assertEqual(utils.row_to_block(row),map[row])

        with self.assertRaises(Exception) as context:
            utils.row_to_block(40)
            self.assertTrue('Too many rows in array' in context.exception)

    def test_all_subsets(self):
        raw_spot_collections = ["2018-01-24_E14_X31",
                                "2018-01-24_E15_X31",
                                "2018-01-24_N21_Pan",
                                "2018-01-24_N22_Cal",
                                "2018-01-24_N23_X31",
                                ]
        ss = itertools.combinations(raw_spot_collections,2)
        print(list(ss))




if __name__ == '__main__':
    unittest.main()
