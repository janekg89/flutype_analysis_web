import copy
from spot2intensity import Point,Rectangle, Grid, Collection

STUDY = "2018-02-21_microarray"
MEASUREMENTS = ["2018-02-21_Bris10_37",
                "2018-02-21_Bris59_26",
                "2018-02-21_Cal_31",
                "2018-02-21_HK_22",
                "2018-02-21_HK_23",
                "2018-02-21_HK_24",
                "2018-02-21_HK_25",
                "2018-02-21_X31_43",
                ]
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

COLLECTIONS["2018-02-21_Bris10_37"]= Collection(grids=grids,
                                          name="2018-02-21_Bris10_37",
                                          jpg_path=FPATH_DICT["2018-02-21_Bris10_37"],
                                          study=STUDY,
                                          pep_path=PEP_PATH,)

#############################