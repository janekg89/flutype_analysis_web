import math
import numpy as np


class Point(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __neg__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __repr__(self):
        return "Point({},{})".format(self.x, self.y)

    def __str__(self):
        return "({},{})".format(self.x, self.y)

    def __abs__(self):
        return Point(abs(self.x), abs(self.y))



    @property
    def as_tuple(self):
        return (self.x, self.y)

    @property
    def abs_distance(self):
        return math.sqrt(abs(self.x) ** 2 + abs(self.y) ** 2)


class Rectangle(object):

    def __init__(self, corner1, corner2, corner3, corner4):
        self.corner1 = corner1
        self.corner2 = corner2
        self.corner3 = corner3
        self.corner4 = corner4

    def corner(self,num):
        map_dict = {1:self.corner1,2:self.corner2,3:self.corner3, 4:self.corner4}
        return map_dict[num]

    @property
    def size(self):
        return self.corner3 - self.corner1

    @property
    def vertical(self):
        return self.corner2 - self.corner1

    @property
    def horizontal(self):
        return self.corner4 - self.corner1

    @property
    def corners(self):
        return [self.corner1.as_tuple,
                self.corner2.as_tuple,
                self.corner3.as_tuple,
                self.corner4.as_tuple,
                ]

    @property
    def corners_as_points(self):
        return [self.corner1,
                self.corner2,
                self.corner3,
                self.corner4,
                ]


    @property
    def rotational_angle(self):
        return math.tan(self.vertical.x / self.vertical.y)


class Grid(object):

    def __init__(self, rectangle, shape):
        self.rectangle = rectangle
        self.shape = shape


    def __str__(self):
        return str(self.rectangle.corner1)

    @property
    def horizontal_points(self, with_off_set=False):
        x = np.linspace(0, self.rectangle.horizontal.x, self.shape.x)
        y = np.linspace(0, self.rectangle.horizontal.y, self.shape.x)
        return list(zip(x, y))

    @property
    def vertical_points(self, with_off_set=False):
        x = np.linspace(0, self.rectangle.vertical.x, self.shape.y)
        y = np.linspace(0, self.rectangle.vertical.y, self.shape.y)
        return list(zip(x, y))

    @property
    def abs_vertical_spacing(self):
        return Point(self.vertical_points[1][0], self.vertical_points[1][1]).abs_distance

    @property
    def abs_horizontal_spacing(self):
        return Point(self.horizontal_points[1][0], self.horizontal_points[1][1]).abs_distance

    @property
    def points(self):

        horizontal_points = [(x + self.rectangle.corner1.x, y + self.rectangle.corner1.y) for x, y in
                             self.horizontal_points]
        xy = []
        for vertical_points_x, vertical_points_y in self.vertical_points:
            xy += [(x + vertical_points_x, y + vertical_points_y) for x, y in horizontal_points]

        return xy

    def add_row(self, where="top"):

        shift = Point(self.vertical_points[1][0], self.vertical_points[1][1])
        self.shape.y += 1

        if where == "top":
            self.rectangle.corner1 = self.rectangle.corner1 - shift
            self.rectangle.corner4 = self.rectangle.corner4 - shift

        if where == "bottom":
            self.rectangle.corner2 = self.rectangle.corner2 + shift
            self.rectangle.corner3 = self.rectangle.corner3 + shift