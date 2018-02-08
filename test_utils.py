import utils
import unittest


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


if __name__ == '__main__':
    unittest.main()
