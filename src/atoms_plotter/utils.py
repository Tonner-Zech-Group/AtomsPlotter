import ase
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np
from ase.data import covalent_radii,colors
from ase.neighborlist import NeighborList
from ase.build import sort
import json
import matplotlib as mpl
from sys import argv
import matplotlib.patheffects as mpe
import math
from ase.build import make_supercell
from matplotlib import collections as mc
from collections import deque
from itertools import islice

def sliding_window(iterable, n):
      """
      sliding_window('ABCDEFG', 4) -> ABCD BCDE CDEF DEFG

      recipe from python docs
      """
      it = iter(iterable)
      window = deque(islice(it, n), maxlen=n)
      if len(window) == n:
          yield tuple(window)
      for x in it:
          window.append(x)
          yield tuple(window)
def color_gradient(x, y, c1=(1,1,1), c2=(0,0,0),zorder=0,linewidth=1,capstyle='projecting',outline=None,NUM=3): #draw gradient bonds
    """
    Creates a line collection with a gradient from colors c1 to c2,
    from data x and y.
    """
    X=np.linspace(x[0],x[1],num=NUM)
    Y=np.linspace(y[0],y[1],num=NUM)
    colors=np.linspace(c1,c2,NUM-1)
    if outline == None:
        return mc.LineCollection(sliding_window(zip(X, Y), 2),
              colors=colors,
              zorder=zorder,
              linewidth=linewidth,
              capstyle=capstyle)
    else:
          return mc.LineCollection(sliding_window(zip(X, Y), 2),
              colors=colors, 
              zorder=zorder,
              linewidth=linewidth,
              capstyle=capstyle,path_effects=[outline])#self.capstyle)#, c1, c2))

def is_inside_cell(pos, cell):
    # Get the atom's position
    atom_pos = pos

    # Invert the cell matrix to get the transformation matrix
    inv_cell = np.linalg.inv(cell)

    # Calculate fractional coordinates of the atom
    fractional_coords = np.dot(atom_pos, inv_cell)

    # Check if the fractional coordinates are within [0, 1] in all dimensions
    return all(0 <= coord <= 1 for coord in fractional_coords)


