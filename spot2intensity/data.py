import copy
from spot2intensity import Point,Rectangle, Grid, Collection


STUDY = "2018-01-24_microarray"
MEASUREMENTS = ["180124_N21_Pan", "180124_E14_X31", "180124_N22_Cal", "180124_N23_X31", "180124_E15_X31"]
DIRECTORY = ["data/{}/{}".format(STUDY,measurement) for measurement in MEASUREMENTS]
IMAGE_PATHS = ["180124_N21_Pan.gif",
                '180124_E14_X31_Allantois_100_600_635_635_2_635.gif',
                '180124_N22_Cal_100_600_635_635_635.gif',
                '180124_N23_X31_100_600_635_635_635_635.gif',
                '180124_E15_X31_100_600_635_635.gif',
               ]

IMAGE_PATHS_JPG  = [im_path.replace('.gif', '.jpg') for im_path in IMAGE_PATHS]
FPATH_DICT = dict(zip(DIRECTORY, IMAGE_PATHS_JPG))
FPATH_DICT = {key: 'data/' + key + '/' + value for key, value in FPATH_DICT.items()}
BLOCK_SHAPE = Point(10,10)
COLLECTIONS = {}
PEP_PATH = "base_pep.gal"

########################################################################


rec_E14_1 = Rectangle(Point(296, 244),
                      Point(314, 1746),
                      Point(1903, 1766),
                      Point(1912, 207)
                      )


rec_E14_2 = Rectangle(Point(303, 2791),
                      Point(298, 4303),
                      Point(1922, 4293),
                      Point(1908, 2786),
                      )

rec_E14_3 = Rectangle(Point(302, 5164),
                      Point(297, 6680),
                      Point(1923, 6701),
                      Point(1910, 5158),
                     )
recs = [rec_E14_1,rec_E14_2,rec_E14_3]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["180124_E14_X31"]= Collection(grids=grids,
                                          name="180124_E14_X31",
                                          jpg_path=FPATH_DICT["180124_E14_X31"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)

#############################

rec_N21_1 = Rectangle(Point(326, 329),
                        Point(321, 1845),
                        Point(1938, 1847),
                        Point(1947, 325),
                      )

rec_N21_2 = Rectangle(Point(315, 2881),
                     Point(321, 4410),
                     Point(1939, 4413),
                     Point(1949, 2875),
                      )

rec_N21_3 = Rectangle(Point(325, 5269),
                     Point(326, 6784),
                     Point(1943, 6777),
                     Point(1952, 5259),
                      )

recs = [rec_N21_1,rec_N21_2,rec_N21_3]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["180124_N21_Pan"]= Collection(grids=grids,
                                          name="180124_N21_Pan",
                                          jpg_path=FPATH_DICT["180124_N21_Pan"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)
##############################
rec_N22_1 = Rectangle(Point(310, 314),
            Point(323, 1846),
            Point(1928, 1831),
            Point(1931, 304),
                      )

rec_N22_2 = Rectangle(Point(321, 2873),
            Point(337, 4399),
            Point(1951, 4382),
            Point(1926, 2869),
                      )
rec_N22_3 = Rectangle(Point(311, 5254),
            Point(349, 6768),
            Point(1944, 6766),
            Point(1959, 5239),
                      )

recs = [rec_N22_1,rec_N22_2,rec_N22_3]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["180124_N22_Cal"]= Collection(grids=grids,
                                          name="180124_N22_Cal",
                                          jpg_path=FPATH_DICT["180124_N22_Cal"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)
#################################
rec_N23_1 = Rectangle(Point(336, 321),
            Point(354, 1831),
            Point(1959, 1853),
            Point(1968, 312),
                      )

rec_N23_2 = Rectangle(Point(328, 2858),
            Point(332, 4385),
            Point(1946, 4403),
            Point(1960, 2873),
                      )

rec_N23_3 = Rectangle(Point(321, 5251),
            Point(310, 6763),
            Point(1929, 6789),
            Point(1941, 5255),
                      )

recs = [rec_N23_1,rec_N23_2,rec_N23_3]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["180124_N23_X31"]= Collection(grids=grids,
                                          name="180124_N23_X31",
                                          jpg_path=FPATH_DICT["180124_N23_X31"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)
#######################
rec_E15_1 = Rectangle(Point(292, 250),
            Point(329, 1760),
            Point(1939, 1747),
            Point(1932, 267),
                      )

rec_E15_2 = Rectangle(Point(297, 2810),
            Point(326, 4315),
            Point(1939, 4309),
            Point(1917, 2795),
                      )

rec_E15_3 = Rectangle(Point(297, 5188),
            Point(317, 6692),
            Point(1911, 6707),
            Point(1925, 5193),
                      )

recs = [rec_E15_1,rec_E15_2,rec_E15_3]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["180124_E15_X31"]= Collection(grids=grids,
                                          name="180124_E15_X31",
                                          jpg_path=FPATH_DICT["180124_E15_X31"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)

