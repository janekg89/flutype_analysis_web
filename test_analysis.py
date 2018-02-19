import analysis as a
from  utils import all_same
import unittest
import sys, os
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_webapp')
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_analysis')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "flutype_webapp.settings")
import django
django.setup()
from flutype.models import Spot
studies = ["2018-01-24_microarray"]


class AnalysisTestCase(unittest.TestCase):

    def setUp(self):
        self.studies = ["2018-01-24_microarray"]
        self.spots = Spot.objects.filter(raw_spot__raw_spot_collection__studies__sid__in=studies)
        self.data = a.Data(self.spots, mean_on=None)

    def test_add_replica_row(self):
        pass

    def test_create_con_diff_features(self):
        spots = self.data._reformat(self.spots)
        spots = self.data._add_replica_row(spots)
        spots = self.data._append_con_diff_features(spots)

        constructed_spots = spots[spots["Constructed Feature"]]
        self.assertTrue(len(constructed_spots.index) > 0)
        self.assertTrue(constructed_spots["Ligand Batch"].count() > 0)
        self.assertTrue(constructed_spots["Ligand Batch Concentration"].count() == 0)


    def test_normalize(self):
        spots = self.data._reformat(self.spots)
        spots = self.data._add_replica_row(spots)
        spots = self.data._append_con_diff_features(spots)
        spots = self.data._append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection")
        normalized_spots = spots[spots["Normalized Feature"]]
        self.assertTrue(len(normalized_spots) > 0)

    def test_pivot_table(self):
        spots = self.data._reformat(self.spots)
        spots = self.data._add_replica_row(spots)
        d_pivot = self.data._pivot_table(spots)
        spots = self.data._append_con_diff_features(spots)
        d_pivot1 = self.data._pivot_table(spots)
        spots = self.data._append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection")
        d_pivot2 = self.data._pivot_table(spots)

        self.assertTrue(all_same([d_pivot.shape[0],d_pivot1.shape[0],d_pivot2.shape[0]]))
        self.assertTrue(d_pivot.shape[1] < d_pivot1.shape[1] < d_pivot2.shape[1])






















if __name__ == '__main__':

    unittest.main()
