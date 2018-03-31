import tkinter as tk
from PIL import ImageTk, Image
from .model import Point, Grid, Rectangle
import numpy as np

class FindGrid(tk.Frame):
    def __init__(self, parent,collection,  *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)

        self.parent = parent
        self.event2canvas = lambda e, c: (c.canvasx(e.x), c.canvasy(e.y))

        self.frame = tk.Frame(parent, bd=2, relief=tk.SUNKEN)
        self.frame.grid_rowconfigure(0, weight=1)
        self.frame.grid_columnconfigure(0, weight=1)

        xscroll = tk.Scrollbar(self.frame, orient=tk.HORIZONTAL)
        xscroll.grid(row=1, column=0, sticky=tk.E + tk.W)
        yscroll = tk.Scrollbar(self.frame)
        yscroll.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.canvas = tk.Canvas(self.frame, bd=0, xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        self.image = self._reformat_image(collection)

        self.canvas.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.canvas.background = self._reformat_image(collection)
        self.canvas.create_image(0, 0, image=self.image, anchor="nw")

        xscroll.config(command=self.canvas.xview)
        yscroll.config(command=self.canvas.yview)
        self.frame.pack(fill=tk.BOTH, expand=1)



        self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))

        self.grids = collection.grids
        self.shape = collection.grids[0].shape
        self.corners = collection.grids[0].rectangle.corners_as_points
        self.collection = collection


        # mouseclick event
        self.canvas.bind_all('<KeyPress-q>', self.left_top_corner)
        self.canvas.bind_all('<KeyPress-w>', self.right_top_corner)
        self.canvas.bind_all('<KeyPress-a>', self.left_bottom_corner)
        self.canvas.bind_all('<KeyPress-s>', self.right_bottom_corner)
        self.canvas.bind_all('<KeyPress-g>', self.show_grid)
        self.canvas.bind_all('<KeyPress-y>', self.save)
        self.canvas.bind_all('<KeyPress-r>', self.reload)


    @staticmethod
    def _reformat_image(collection):
        try:
            im = collection.tifs_a[1]
        except:
            im = collection.tifs_a[0]

        im_array = np.asarray(im)
        im_array = im_array * (255.0 / im_array.max())
        return ImageTk.PhotoImage(Image.fromarray(im_array))

    def point2grid(self, i, event):
        cx, cy = self.event2canvas(event, self.canvas)
        self.corners[i] = Point(cx, cy)
        print ("Point(%d, %d)" % (cx, cy))
        return Point(cx, cy)

    def left_top_corner(self,event):
        self.point2grid(0,event)
        self.canvas.create_image(0, 0, image=self._reformat_image(self.collection), anchor="nw")

    def left_bottom_corner(self,event):
        self.point2grid(1,event)

    def right_bottom_corner(self,event):
        self.point2grid(2, event)

    def right_top_corner(self, event):
        self.point2grid(3,event)

    def show_grid(self,event):
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, image=self.image, anchor="nw")
        this_grid = Grid(Rectangle(*self.corners),self.shape)
        for x, y in this_grid.points:
            delta_x = int(this_grid.abs_horizontal_spacing)
            delta_y = int(this_grid.abs_vertical_spacing)
            x = int(x)
            y = int(y)
            self.canvas.create_rectangle(x-.5*delta_x,
                                    y-.5*delta_y,
                                    x+.5*delta_x,
                                    y+.5*delta_y,
                                    outline="green")

    def save(self,event):
        pass

    def reload(self,event):
        pass



