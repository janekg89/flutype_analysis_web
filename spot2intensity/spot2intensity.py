import math
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
from skimage.transform import hough_circle, hough_circle_peaks
from skimage import feature
import skimage.io
from PIL import Image, ImageSequence
import attr
import pickle
from Tkinter import *
from PIL import ImageTk, Image


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

@attr.s
class Spots(object):
    df = attr.ib()
    b_im = attr.ib(default=None)
    a_im = attr.ib(default=None)


    def plot_grid(self, on="a_im", **kwargs):
        image = getattr(self,on)
        if not kwargs.get("ax"):
            fig, ax = plt.subplots(1, **kwargs)
        else:
            ax = kwargs.get("ax")

        ax.imshow(np.asarray(image), cmap="gray")
        for i, spot in self.df.iterrows():
            circ = create_circle_patches(spot["circles"].center, spot["circles"].radius)
            rec = patches.Rectangle((spot["squares"].get_xy()),
                                    spot["squares"].get_width(),
                                    spot["squares"].get_height(),
                                    fill=False,  # remove background
                                    linewidth=1,
                                    edgecolor='g',
                                    )
            ax.add_patch(circ)
            ax.add_patch(rec)

    @staticmethod
    def load_pickel(collection):
        directory = "data/{}/{}/".format(collection.study, collection.name)
        with open(directory + "spots_class", 'r') as f:
            spots = pickle.load(f)
        return spots

    def select_by_circlequal(self,low_tresh,inplace=True):
        selected = self.df.loc[self.df["circle_qual"] > low_tresh * self.df["circle_qual"].max()]
        if not inplace:
            return selected
        self.df = selected


    def add_virus(self,virus):
        self.df["Virus"] = virus

    def add_c_name(self, c_name):
        self.df["Collection"] = c_name


class Collection(object):
    def __init__(self, grids,  name, study, jpg_path=None,tif_a_path=None,tif_b_path=None,pep_path=None,grid_a_path=None,grid_b_path=None):
        self.study = study
        self.name = name
        self.grids = grids
        self.png_path = jpg_path
        self.image = skimage.io.imread(jpg_path, 0) if jpg_path is not None else None
        self.tif_a_path = tif_a_path
        self.tif_b_path = tif_b_path
        self.pep_path = pep_path
        self.gal = pd.read_csv(pep_path, sep='\t', index_col="ID") if pep_path is not None else None
        self.r,self.g,self.b = self.color_splitted_image() if jpg_path is not None else (None,None,None)
        self.tifs_a = self.read_tif(self.tif_a_path) if tif_a_path is not None else [None,None]
        self.tifs_b = self.read_tif(self.tif_b_path) if tif_b_path is not None else [None,None]
        self.grid_a_path = grid_a_path
        self.grid_b_path = grid_b_path

    @staticmethod
    def read_tif(tif):
        """
        path - Path to the multipage-tiff file
        n_images - Number of pages in the tiff file
        """
        img = Image.open(tif)
        images = []
        for i, page in enumerate(ImageSequence.Iterator(img)):
            images.append(page.copy())
        return images

    def color_splitted_image(self):
        r, g, b = np.dsplit(self.image, 3)

        r = r.squeeze()
        g = g.squeeze()
        b = b.squeeze()
        return r,g,b

    def create_values(self,of="tif_a"):
        """
        of can be "tif_b",tif_a"
        :param of:
        :return:
        """
        spot_images = []
        intensities2_b = []
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
                  "intensities2_b":intensities2_b,
                  "std_intensities2":std_intensities2,
                  "circle_qual":circle_qual,
                  }

        images_dict = {"tif_a": np.asarray(self.tifs_a[1]), "tif_b": np.asarray(self.tifs_b[1])}
        image = images_dict.get(of)
        if of == "tif_a":
            image_b = images_dict["tif_b"]
        elif of == "tif_b":
            image_b = images_dict["tif_a"]

        print(self.grids)
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

                spot_imag = image[y0:y0 + delta_x, x0:x0 + delta_y]
                spot_imag_b = image_b[y0:y0 + delta_x, x0:x0 + delta_y]

                spot_images.append(spot_imag)
                intensities.append(spot_imag.sum())

                circx, circy, radius, accums = find_circle_coordinates(spot_imag)
                circle_qual.append(accums)
                circ = create_circle_patches((circy + x0, circx + y0), radius)
                circles.append(circ)

                circ_points = np.array([value for (y, x), value in np.ndenumerate(spot_imag) if
                                        contained_in_circle(circx, x, circy, y, radius)])
                circ_points_b = np.array([value for (y, x), value in np.ndenumerate(spot_imag_b) if
                                        contained_in_circle(circx, x, circy, y, radius)])

                std_intensities2.append(circ_points.std())
                intensities2.append(circ_points.mean())
                intensities2_b.append(circ_points_b.mean())


        return values

    def pd_complete_spots(self, **kwargs):

        spots = self.gal
        values = self.create_values(**kwargs)
        for key,value in values.items():
            spots[key] = value
        im_b = np.asarray(self.tifs_b[1])
        im_a = np.asarray(self.tifs_a[1])

        return Spots(df=spots, b_im=im_b, a_im =im_a)

    def create_grid(self,shape,on="tiff_a"):
        tiff_dict = {"tiff_a":self.tifs_a[1],"tiff_b":self.tifs_b[1]}
        assert tiff_dict[on] is not None, " You have no image loaded"
        im_array = np.asarray(tiff_dict[on])
        im_array = im_array * (255.0 / im_array.max())
        img = Image.fromarray(im_array)
        self.canvas_microarray(img,shape)
        #assert all(self.tiff_a[1] is None, self.tiff_b[1]  is None) , "no images loaded"

    def canvas_microarray(self, image, shape):
        event2canvas = lambda e, c: (c.canvasx(e.x), c.canvasy(e.y))
        root = Toplevel()
        img = ImageTk.PhotoImage(image)
        # setting up a tkinter canvas with scrollbars
        frame = Frame(root, bd=2, relief=SUNKEN)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        xscroll = Scrollbar(frame, orient=HORIZONTAL)
        xscroll.grid(row=1, column=0, sticky=E + W)
        yscroll = Scrollbar(frame)
        yscroll.grid(row=0, column=1, sticky=N + S)
        canvas = Canvas(frame, bd=0, xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        canvas.grid(row=0, column=0, sticky=N + S + E + W)
        xscroll.config(command=canvas.xview)
        yscroll.config(command=canvas.yview)
        frame.pack(fill=BOTH, expand=1)
        canvas.create_image(0, 0, image=img, anchor="nw")
        canvas.config(scrollregion=canvas.bbox(ALL))
        # function to be called when mouse is clicked
        grids = [grid.rectangle.corners_as_points for grid in self.grids]
        grid = grids[0]

        def point2grid(i,event):
            cx, cy = event2canvas(event, canvas)
            print ("Point(%d, %d)" % (cx, cy))
            return Point(cx, cy)

            #grid[i] = Point(cx, cy)

        def left_top_corner(event):
            # outputting x and y coords to console
            #point2grid(0,event)
            cx, cy = event2canvas(event, canvas)
            print ("Point(%d, %d)" % (cx, cy))
            #grid[0] = Point(cx, cy)
            return Point(cx, cy)

        def left_bottom_corner(event):
            #point2grid(1,event)
            cx, cy = event2canvas(event, canvas)
            print ("Point(%d, %d)" % (cx, cy))
            grid[1] = Point(cx, cy)

        def right_bottom_corner(event):
            #point2grid(3, event)
            cx, cy = event2canvas(event, canvas)
            print ("Point(%d, %d)" % (cx, cy))
            grid[2] = Point(cx, cy)

        def right_top_corner(event):
            #point2grid(4,event)
            cx, cy = event2canvas(event, canvas)
            print ("Point(%d, %d)" % (cx, cy))
            grid[3] = Point(cx, cy)


        def show_grid(grid,shape):
            this_grid = Grid(Rectangle(*grid),shape)
            canvas.delete("all")
            canvas.create_image(0, 0, image=img, anchor="nw")
            for x, y in this_grid.points:
                delta_x = int(this_grid.abs_horizontal_spacing)
                delta_y = int(this_grid.abs_vertical_spacing)
                x = int(x)
                y = int(y)
                canvas.create_rectangle(x-.5*delta_x,
                                        y-.5*delta_y,
                                        x+.5*delta_x,
                                        y+.5*delta_y,
                                        outline="green" )


        # mouseclick event
        canvas.bind_all('<KeyPress-q>', left_top_corner)
        canvas.bind_all('<KeyPress-a>', left_bottom_corner)
        canvas.bind_all('<KeyPress-w>', right_top_corner)
        canvas.bind_all('<KeyPress-s>', left_bottom_corner)
        canvas.bind_all('<KeyPress-q>', show_grid(grid,shape))

        root.mainloop()

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


def find_circle_coordinates(image):

    pic = copy.deepcopy(image)


    #pic[pic > 120] = 0
    #pic[pic < 1] = 0
    #edges = roberts(pic)
    edges = feature.canny(pic)


    # Detect two radii
    hough_radii = np.arange(10, 60, 2)

    hough_res = hough_circle(edges, hough_radii)

    # Select the most prominent  circle:
    accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, min_xdistance=hough_radii.min(),
                                               min_ydistance=hough_radii.min(),
                                               total_num_peaks=1)
    if len(cx) == 0 :
        radii = [20]
        accums = [0]
        x,y = image.shape

        cy = [y/2]
        cx = [x/2]


    return cy[0] , cx[0], radii[0], accums[0]