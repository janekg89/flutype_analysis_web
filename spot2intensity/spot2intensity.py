import math
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
from skimage.transform import hough_circle, hough_circle_peaks
from skimage import feature
from skimage.filters import roberts
import cv2
import skimage.io



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


class Collection(object):

    def __init__(self, grids, jpg_path,pep_path, name, study):
        self.study = study
        self.name = name
        self.grids = grids
        self.png_path = jpg_path
        self.image = skimage.io.imread(jpg_path, 0)
        self.pep_path = pep_path
        self.gal = pd.read_csv(pep_path, sep='\t', index_col="ID")
        self.r,self.g,self.b = self.color_splitted_image()

    def color_splitted_image(self):
        r, g, b = np.dsplit(self.image, 3)

        r = r.squeeze()
        g = g.squeeze()
        b = b.squeeze()
        return r,g,b

    def create_values(self):
        spot_images = []
        intensities = []
        intensities2 = []
        circles = []
        squares = []
        circle_qual = []
        std_intensities2 = []

        values = {"spot_images":spot_images,
                  "squares": squares,
                  "circles": circles,
                  "intensities":intensities,
                  "intensities2":intensities2,
                  "std_intensities2":std_intensities2,
                  "circle_qual":circle_qual,
                  }

        for grid in self.grids:
            print("Claculating grid starting at:{}".format(grid))

            #grid.add_row(where="top")
            #grid.add_row(where="bottom")

            for x, y in grid.points:
                delta_x = int(grid.abs_horizontal_spacing)
                delta_y = int(grid.abs_vertical_spacing)
                x = int(x)
                y = int(y)
                rec = create_patches((x, y), delta_x, delta_y)

                squares.append(rec)
                x0, y0 = rec.xy
                x0 = int(x0)
                y0 = int(y0)
                if x0 < 0:
                    x0 = 0
                if y0 < 0:
                    y0 = 0

                spot_imag = self.r[y0:y0 + delta_x, x0:x0 + delta_y]
                spot_images.append(spot_imag)
                intensities.append(spot_imag.sum())

                circx, circy, radius, accums = find_circle_coordinates(spot_imag, 0, 0)
                circle_qual.append(accums)
                circ = create_circle_patches((circy + x0, circx + y0), radius)
                circles.append(circ)

                circ_points = np.array([value for (y, x), value in np.ndenumerate(spot_imag) if
                                        contained_in_circle(circx, x, circy, y, radius)])

                std_intensities2.append(circ_points.std())
                intensities2.append(circ_points.mean())

        return values

    def pd_complete_spots(self):

        spots = self.gal
        values = self.create_values()
        for key,value in values.iteritems():
            spots[key] = value
        return spots



def contained_in_circle(x0,x,y0,y,r):
    return r**2 >= (x0-x)**2 + (y0-y)**2

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
        edgecolor='b',
    )

    return rec


def rectangle_reshape(x_0, y_0, x_spacing, y_spacing):
    x_min = int(x_0 - 0.5 * x_spacing)
    x_max = int(x_0 + 0.5 * x_spacing)
    y_min = int(y_0 - 0.5 * y_spacing)
    y_max = int(y_0 + 0.5 * y_spacing)
    return x_min, x_max, y_min, y_max


def find_circle_coordinates(image, x0, y0):

    pic = copy.deepcopy(image)


    #pic[pic > 120] = 0
    #pic[pic < 1] = 0
    edges = roberts(image)
    #edges = feature.canny(pic)

    # Detect two radii
    hough_radii = np.arange(10, 60, 2)

    hough_res = hough_circle(edges, hough_radii)

    # Select the most prominent  circle:
    accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, min_xdistance=hough_radii.min(),
                                               min_ydistance=hough_radii.min(),
                                               total_num_peaks=1)
    return cy[0] + y0, cx[0] + x0, radii[0], accums[0]