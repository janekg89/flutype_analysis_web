import os
import copy
import attr
import pickle
import datetime

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from skimage.transform import hough_circle, hough_circle_peaks
from skimage import feature, exposure, filters
import skimage.io
from PIL import Image, ImageSequence

from .ui_grid import FindGrid
from .model import Point, Rectangle, Grid
from utils import ensure_dir


try:
    from six.moves import tkinter as Tk
    # print('Python 3')

except ImportError :
    from six.moves import Tkinter as Tk
    # print 'Python 2.7'

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

        ax.imshow(exposure.equalize_hist(np.asarray(image)), cmap="gray")
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
        directory = os.path.join(collection.base_path, collection.name)
        #directory = ( collection.name)
        with open(os.path.join(directory, "spots_class.pickel"), 'rb') as f:
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

    def raw_intensies_pivot(self):
        return self.df.pivot(index="Row",columns="Column",values="intensities")


    def quant1_intensies_pivot(self):
        return self.df.pivot(index="Row",columns="Column",values="intensities2")

    def before_intensies_pivot(self):
        return self.df.pivot(index="Row",columns="Column",values="intensities2_b")

    def std_intensies2_pivot(self):
        return self.df.pivot(index="Row",columns="Column",values="std_intensities2")

    def circle_qual_intensies2_pivot(self):
        return self.df.pivot(index="Row",columns="Column",values="circle_qual")



class Collection(object):
    def __init__(self, grids,  name, study,base, jpg_path=None,tif_a_path=None,tif_b_path=None,pep_path=None,grid_a_path=None,grid_b_path=None, virus =""):
        self.base_path = base
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
        self.collection_meta = self.c_meta_dict()
        self.collection_steps = self.c_steps_pd()
        self.virus = virus
        self.gal_vir = self.gal_vir()
        self.raw_meta = self.r_meta_dict_raw()
        self.quant1_meta = self.r_meta_dict_quant1()
        self.before_meta = self.r_meta_dict_before()


    def gal_vir(self):
        gal_vir = copy.deepcopy(self.gal)
        gal_vir["Name"] = self.virus
        return gal_vir

    def resize_tif_b(self):
        #im = Image.fromarray(np.asarray(self.tifs_b[-1]).astype(np.uint8))
        #im_thumb = im.thumbnail(self.tifs_a[-1].size, Image.ANTIALIAS)
        #output = self.base_path+"tif_b_new.tif"
        #im.save(output)
        return "hi"

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

        images_dict = {"tif_a": np.asarray(self.tifs_a[-1]),
                       "tif_b": np.asarray(self.tifs_b[-1])}
        image = images_dict.get(of)
        if of == "tif_a":
            image_b = images_dict["tif_b"]
        elif of == "tif_b":
            image_b = images_dict["tif_a"]

        equilized_image = exposure.equalize_hist(image)

        print(self.grids)
        for grid in self.grids:
            print("Claculating grid starting at:{}".format(grid))

            #grid.add_row(where="top")
            #grid.add_row(where="bottom")
            delta_x = int(grid.abs_horizontal_spacing)
            delta_y = int(grid.abs_vertical_spacing)

            for x, y in grid.points:

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
                spot_imag_eq = equilized_image[y0:y0 + delta_x, x0:x0 + delta_y]


                spot_imag_b = image_b[y0:y0 + delta_x, x0:x0 + delta_y]

                spot_images.append(spot_imag)
                intensities.append(spot_imag.sum())

                circy, circx, radius, accums = find_circle_coordinates(spot_imag_eq)
                circle_qual.append(accums)
                circ = create_circle_patches((circx + x0, circy + y0), radius)
                circles.append(circ)

                circ_points = np.array([value for (y, x), value in np.ndenumerate(spot_imag) if
                                        contained_in_circle(circy, x, circx, y, radius)])
                circ_points_b = np.array([value for (y, x), value in np.ndenumerate(spot_imag_b) if
                                        contained_in_circle(circy, x, circx, y, radius)])

                std_intensities2.append(circ_points.std())
                intensities2.append(circ_points.mean())
                intensities2_b.append(circ_points_b.mean())


        return values

    def pd_complete_spots(self, **kwargs):

        spots = self.gal
        values = self.create_values(**kwargs)
        for key,value in values.items():
            spots[key] = value
        im_b = np.asarray(self.tifs_b[-1])
        im_a = np.asarray(self.tifs_a[-1])

        return Spots(df=spots, b_im=im_b, a_im =im_a)

    def dump_pickel(self):
        print("<{}> is beeing pickeled".format(self.name))
        spots = self.pd_complete_spots()
        with open(os.path.join(self.base_path,self.name, "spots_class.pickel"), 'wb') as f:
                pickle.dump(spots,f)


    def create_grid(self,on="tiff_a"):
        tiff_dict = {"tiff_a":self.tifs_a[-1],"tiff_b":self.tifs_b[-1]}
        assert tiff_dict[on] is not None, "You have no image loaded"
        root = Tk.Toplevel()
        app = FindGrid(root,collection=self)
        root.mainloop()
        print(app.corners)

    def export_measurement_folder(self,where):
        print("Starting with Collection <{}>".format(self.name))
        assert os.path.exists(where), "<{}> destination does not exist".format(where)
        #creates all nessesary folders
        collection_folder = os.path.join(where,self.name)
        results_folder = os.path.join(collection_folder,"results")

        raw_folder = os.path.join(results_folder,"raw/")
        quant1_folder = os.path.join(results_folder,"quant1/")
        before_folder = os.path.join(results_folder,"before/")

        ensure_dir(raw_folder)
        ensure_dir(quant1_folder)
        ensure_dir(before_folder)

        # main folder
        pd.DataFrame.from_dict(data=self.collection_meta, orient='index').to_csv(os.path.join(collection_folder,"meta.tsv"), sep=str('\t'), header=False)

        self.c_steps_pd().to_csv(os.path.join(collection_folder,"steps.tsv"), sep=str('\t'), encoding='utf-8')
        self.gal.to_csv(os.path.join(collection_folder,"lig_fix.txt"), sep=str('\t'), index="True", encoding='utf-8')
        self.gal_vir.to_csv(os.path.join(collection_folder,"lig_mob.txt"), sep=str('\t'), index="True", encoding='utf-8')

        spots = Spots.load_pickel(self)

        #Image.fromarray(spots.b_im).save(os.path.join(collection_folder,"b_im.tiff"),"TIFF")
        #Image.fromarray(spots.a_im).save(os.path.join(collection_folder,"a_im.tiff"),"TIFF")
        Image.fromarray(spots.b_im).point(lambda i:i*(1./256)).convert('L').save(os.path.join(collection_folder, "b_im.png"))
        Image.fromarray(spots.a_im).point(lambda i:i*(1./256)).convert('L').save(os.path.join(collection_folder, "a_im.png"))

        #except:
        #    spots = self.pd_complete_spots()

        #results before
        pd.DataFrame.from_dict(data=self.before_meta, orient='index').to_csv(os.path.join(before_folder,"meta.tsv"), sep=str('\t'), header=False)
        spots.before_intensies_pivot().to_csv(os.path.join(before_folder,"intensity.tsv"), sep=str('\t'))

        # results raw
        pd.DataFrame.from_dict(data=self.raw_meta, orient='index').to_csv(os.path.join(raw_folder,"meta.tsv"), sep=str('\t'), header=False)
        spots.raw_intensies_pivot().to_csv(os.path.join(raw_folder,"intensity.tsv"), sep=str('\t'))

        # results quant1
        pd.DataFrame.from_dict(data=self.quant1_meta, orient='index').to_csv(os.path.join(quant1_folder,"meta.tsv"), sep=str('\t'), header=False)
        spots.quant1_intensies_pivot().to_csv(os.path.join(quant1_folder,"intensity.tsv"), sep=str('\t'))
        spots.circle_qual_intensies2_pivot().to_csv(os.path.join(quant1_folder,"circle_quality.tsv"), sep=str('\t'))
        spots.std_intensies2_pivot().to_csv(os.path.join(quant1_folder,"std.tsv"), sep=str('\t'))


    def r_meta_dict_raw(self):
        d = {"comment":
        "spotdetection with spot2intensity.py from {} .".format(datetime.date.today()),
        "processing_type": "Fluorescence reader of microwell plate (sum over square on grid )",
        "intensity_file": "intensity.tsv"}
        return d

    def r_meta_dict_quant1(self):
        d = self.r_meta_dict_raw()
        d["processing_type"] = "Fluorescence reader of microwell plate (mean on best fitfing circle of spot)"
        d["circle_quality"] = "circle_quality.tsv"
        d["std"] = "std.tsv"
        return d

    def r_meta_dict_before(self):
        d = self.r_meta_dict_raw()
        d["comment"] = "spotdetection with spot2intensity.py from {} after spotting.".format(datetime.date.today())
        return d


    def c_steps_pd(self):
        d = [
            ["S03", "0", "", "", "", "", ""],
            ["SC001", "1", "", "", "", "b_im.png", ""],
            ["W004", "2", "", "", "", "", ""],
            ["I01", "3", "", "", "", "", ""],
            ["W004", "4", "", "", "", "", ""],
            ["SC001", "5", "", "", "", "a_im.png", ""],
        ]

        return pd.DataFrame(d, columns=["step", "index", "user", "start", "comment", "image", "intensities"])


    @staticmethod
    def c_meta_dict(func="3D-Epoxy",mt="microarray", ma = "PolyAn", com = ""):
            meta_dict = {"functionalization":func,
                         "measurement_type":mt,
                         "manufacturer":ma,
                         "comment":com}
            return meta_dict



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
        linewidth=2,
        edgecolor='r',
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
    #img_seg_eq = segmentation.felzenszwalb(pic, scale=0.1, sigma=1, min_size=100)
    edges = feature.canny(filters.gaussian(pic, sigma=4))
    #edges = feature.canny(pic)
    # Detect two radii
    hough_radii = np.arange(6, 40, 2)
    hough_res = hough_circle(edges[5:-5,5:-5], hough_radii,full_output=False)
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