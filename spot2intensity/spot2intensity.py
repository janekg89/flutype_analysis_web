import math
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




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
    def rotational_angle(self):
        return math.tan(self.vertical.x / self.vertical.y)


class Grid(object):

    def __init__(self, rectangle, shape):
        self.rectangle = rectangle
        self.shape = shape

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


def plot_patches(x, y, x_spacing, y_spacing):
    fig2 = plt.figure(figsize=(30, 10))
    ax2 = fig2.add_subplot(111, aspect='equal')
    pts = np.vstack([x, y]).reshape(2, -1).T
    for p in [
        patches.Rectangle(
            (patchcenter[1] - 0.5 * x_spacing, patchcenter[0] - 0.5 * y_spacing),
            x_spacing,
            y_spacing,
            fill=False,  # remove background
            linewidth=1,
            edgecolor='r'
        ) for patchcenter in pts
    ]:
        ax2.add_patch(p)

    return fig2


def create_patches(pt, x_spacing, y_spacing):
    rec = patches.Rectangle(
        (pt[0] - 0.5 * x_spacing, pt[1] - 0.5 * y_spacing),
        x_spacing,
        y_spacing,
        fill=False,  # remove background
        linewidth=1,
        edgecolor='r'
    )

    return rec

def create_circle_patches(center, radius):
    rec = patches.Circle(
        xy=center,
        radius=radius,
        fill=False,  # remove background
        linewidth=1,
        edgecolor='r',
    )

    return rec

def rectangle_reshape(x_0, y_0, x_spacing, y_spacing):
    x_min = int(x_0 - 0.5 * x_spacing)
    x_max = int(x_0 + 0.5 * x_spacing)
    y_min = int(y_0 - 0.5 * y_spacing)
    y_max = int(y_0 + 0.5 * y_spacing)
    return x_min, x_max, y_min, y_max