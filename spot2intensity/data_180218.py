import copy
from spot2intensity import Point,Rectangle, Grid, Collection

STUDY = "2018-02-18_microarray"


BLOCK_SHAPE = Point(15,57)
COLLECTIONS = {}
PEP_PATH = "data/2018-02-18_microarray/base_pep.gal"
########################################################################
fpath = "data/2018-02-18_microarray/slide_37/afterIncubation/37_80_280_rot_635.jpg"
rec = Rectangle(Point(43, 98),
                      Point(48, 3441),
                      Point(988, 3435),
                      Point(987, 106)
                      )



recs = [rec]
grids = [Grid(rectangle = rec, shape = copy.deepcopy(BLOCK_SHAPE)) for rec in recs]

COLLECTIONS["slide_37"]= Collection(grids=grids,
                                    name="slide_37",
                                    jpg_path=fpath,
                                    study=STUDY,
                                    pep_path=PEP_PATH,
                                    )

#############################