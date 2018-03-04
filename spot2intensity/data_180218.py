import copy
from spot2intensity import Point,Rectangle, Grid, Collection

BLOCK_SHAPE = Point(15, 57)
COLLECTIONS = {}

def c_data(numb):
    data_dict = {}
    data_dict["name"] = "slide_{}".format(numb)
    data_dict["jpg_path"] = "data/2018-02-18_microarray/slide_{}/afterIncubation/{}_80_600_rot_635.jpg".format(numb,numb)
    data_dict["tif_a_path"] = "data/2018-02-18_microarray/slide_{}/afterIncubation/{}_80_600_rot.tif".format(numb,numb)
    data_dict["jpg_path"] = "data/2018-02-18_microarray/slide_{}/afterIncubation/{}_80_600_rot_635.jpg".format(numb,
                                                                                                               numb)
    data_dict["tif_b_path"] = "data/2018-02-18_microarray/slide_{}/slide_{} 600 80 alle Farben.tif".format(numb,numb)
    data_dict["grid_a_path"] = "data/2018-02-18_microarray/slide_{}/im_grid_after_600_80_635.jpg".format(numb, numb)
    data_dict["grid_b_path"] = "data/2018-02-18_microarray/slide_{}/im_grid_before_600_80_635.jpg".format(numb,numb)
    data_dict["study"] =  "2018-02-18_microarray"
    data_dict["pep_path"] = "data/2018-02-18_microarray/base_pep.gal"
    return data_dict

def c_grids(recs):
    return [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]
def load_data():

    ########################################################################
    recs = [Rectangle(
        Point(46, 96),
        Point(79, 3398),
        Point(996, 3401),
        Point(985, 99),
                    ),]
    grids = c_grids(recs)
    n = 47
    data_47 = c_data(n)
    data_47["tif_a_path"] = None
    data_47["jpg_path"] = None
    data_47["grid_a_path"] = None
    data_47["grid_b_path"] = None

    COLLECTIONS[c_data(n)["name"]] = Collection(grids=grids, **data_47)
    ########################################################################
    recs = [Rectangle(
        Point(54, 151),
        Point(69, 3453),
        Point(994, 3464),
        Point(989, 152),
                    ),]
    grids = c_grids(recs)
    n = 46
    data_46 = c_data(n)
    data_46["tif_a_path"] = None
    data_46["jpg_path"] = None
    data_46["grid_a_path"] = None
    data_46["grid_b_path"] = None
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **data_46)
    ########################################################################
    recs = [Rectangle(
        Point(48, 106),
        Point(65, 3409),
        Point(993, 3414),
        Point(987, 103),
                    ),]
    grids = c_grids(recs)
    n = 45
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(53, 95),
        Point(66, 3401),
        Point(996, 3402),
        Point(999, 100),
                    ),]
    grids = c_grids(recs)
    n = 44
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(55, 84),
        Point(81, 3399),
        Point(1008, 3390),
        Point(988, 36),
                    ),]
    grids = c_grids(recs)
    n = 42
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(51, 89),
        Point(55, 3423),
        Point(992, 3414),
        Point(998, 66),
                          ),]
    grids = c_grids(recs)
    n = 41
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################

    recs = [Rectangle(
        Point(53, 83),
        Point(75, 3382),
        Point(994, 3382),
        Point(989, 82),
                    ),]
    grids = c_grids(recs)
    n = 40
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################

    recs = [Rectangle(
        Point(55, 84),
        Point(81, 3399),
        Point(1008, 3390),
        Point(988, 36),
                    ),]
    grids = c_grids(recs)
    n = 39
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(41, 110),
        Point(87, 3448),
        Point(1013, 3441),
        Point(989, 86),
                    ),]
    grids = c_grids(recs)
    n = 38
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(43, 98),
                    Point(48, 3441),
                    Point(988, 3435),
                    Point(987, 106)
                    ),]
    grids = c_grids(recs)
    n = 37
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(43, 102),
                          Point(46, 3431),
                          Point(986, 3436),
                          Point(991, 108)
                          ),]
    grids = c_grids(recs)
    n=36
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(44, 116),
                          Point(46, 3446),
                          Point(991, 3445),
                          Point(993, 124)
                          ),]
    grids = c_grids(recs)
    n=35
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(35, 98),
                          Point(41, 3425),
                          Point(982, 3430),
                          Point(991, 96)
                          ),]
    grids = c_grids(recs)
    n=34
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(51, 89),
                        Point(55, 3423),
                        Point(992, 3414),
                        Point(998, 66),
                          ),]
    grids = c_grids(recs)
    n=33
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(51, 85),
                        Point(54, 3428),
                        Point(995, 3407),
                        Point(1000, 90),
                          ),]
    grids = c_grids(recs)
    n=32
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(Point(37, 90),
                    Point(58, 3401),
                    Point(999, 3415),
                    Point(979, 84)
                      ,),]
    grids = c_grids(recs)
    n=31
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
                    Point(47, 76),
                    Point(73, 3429),
                    Point(1019, 3433),
                    Point(982, 78),)]
    grids = c_grids(recs)
    n=30
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(39, 152),
        Point(37, 3491),
        Point(984, 3480),
        Point(982, 162),)]
    grids = c_grids(recs)
    n=29
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(39, 155),
        Point(55, 3495),
        Point(988, 3487),
        Point(979, 156),)]
    grids = c_grids(recs)
    n=28
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(42, 144),
        Point(43, 3485),
        Point(988, 3464),
        Point(984, 144),)]
    grids = c_grids(recs)
    n=27
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(35, 148),
        Point(42, 3483),
        Point(973, 3482),
        Point(983, 160),)]
    grids = c_grids(recs)
    n=26
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(35, 152),
        Point(44, 3467),
        Point(977, 3481),
        Point(985, 164),)]
    grids = c_grids(recs)
    n=25
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(43, 144),
        Point(33, 3485),
        Point(957, 3475),
        Point(983, 155),)]
    grids = c_grids(recs)
    n=24
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(39, 150),
        Point(49, 3485),
        Point(980, 3475),
        Point(982, 156),)]
    grids = c_grids(recs)
    n=23
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(39, 134),
        Point(49, 3462),
        Point(968, 3463),
        Point(979, 141),)]
    grids = c_grids(recs)
    n=22
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(55, 132),
        Point(61, 3460),
        Point(997, 3460),
        Point(994, 136),)]
    grids = c_grids(recs)
    n=21
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **c_data(n))
    ########################################################################
    recs = [Rectangle(
        Point(40, 103),
        Point(48, 3469),
        Point(955, 3469),
        Point(975, 116),)]
    grids = c_grids(recs)
    n=20
    data_20 = c_data(n)
    data_20["tif_a_path"] = None
    data_20["jpg_path"] = None
    data_20["grid_a_path"] = None
    data_20["grid_b_path"] = None
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **data_20)
    ########################################################################
    recs = [Rectangle(
        Point(30, 107),
        Point(41, 3473),
        Point(966, 3468),
        Point(973, 115),)]
    grids = c_grids(recs)
    n=19
    data_19 = c_data(n)
    data_19["tif_a_path"] = None
    data_19["jpg_path"] = None
    data_19["grid_a_path"] = None
    data_19["grid_b_path"] = None
    COLLECTIONS[c_data(n)["name"]]= Collection(grids=grids, **data_19)

    return COLLECTIONS