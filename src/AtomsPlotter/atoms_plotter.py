import ase
import matplotlib.pyplot as plt
import numpy as np
from ase.data import covalent_radii, colors
from ase.neighborlist import NeighborList
from ase.build import sort
import matplotlib.patheffects as mpe
import math
from .utils import color_gradient


class atoms_plotter():
    def __init__(self,
                 atoms=ase.Atoms(None),
                 ATOMS=True,
                 view=2,
                 repeat=[1, 1],
                 show_unit_cell=False,
                 bond_radius=1,
                 unit_cell_linestyle='dashed',
                 # I determined this to be a nice size. adjust if necessary.
                 scale=100,
                 show=True,
                 show_bonds=True,
                 use_bondorders=True,
                 double_bond_offset=[0.1, 0.1, 0],
                 triple_bond_offset=[0.07, 0.07, 0],
                 double_bond=1.4,
                 outline_width=1.5,
                 azim=0,
                 elev=0,
                 shift=[0, 0],
                 triple_bond=1.34,
                 colorbonds=True,
                 lewis=False,
                 bondlinewidth=1,
                 name='structure',
                 dimension='2D',
                 constraints=True,
                 bond_gradient=True,
                 draw_outline=True,
                 projectiontype='ortho',
                 format='svg'):

        self.format = format
        self.atoms = atoms
        self.ATOMS = ATOMS
        self.view = view
        self.lewis = lewis
        self.repeat = repeat
        # self.plot_atom_cutoff = plot_atom_cutoff
        # self.plot_atom_colorscaling = plot_atom_colorscaling
        self.show_unit_cell = show_unit_cell
        self.bond_radius = bond_radius
        self.unit_cell_linestyle = unit_cell_linestyle
        # I determined this to be a nice size. adjust if necessary.
        self.scale = scale
        self.show = show
        self.show_bonds = show_bonds
        self.use_bondorders = use_bondorders
        self.bondatoms = []
        self.double_bond_offset = double_bond_offset
        self.triple_bond_offset = triple_bond_offset
        self.outline_width = 1.5
        self.azim = azim
        self.elev = elev
        self.double_bond = double_bond
        self.shift = shift
        self.triple_bond = triple_bond
        self.colorbonds = colorbonds
        self.constraints = constraints
        self.bondlinewidth = bondlinewidth
        self.name = name
        self.dimension = dimension
        self.projectiontype = projectiontype
        self.draw_outline = draw_outline
        self.bond_gradient = bond_gradient
        self.color_dict = {'H': (1.00, 1.00, 1.00),
                           'C': (0.35, 0.35, 0.35),
                           'Si': (0.02, 0.38, 0.67),
                           'Ge': (0.05, 0.45, 0.45),
                           'Ga': (0.33, 0.71, 0.09),
                           'In': (0, 0, 0),
                           'N': (0.00, 0.00, 1.00),
                           'P': (1.00, 0.50, 0.00),
                           'As': (0.75, 0.54, 0.00),
                           'Sb': (0.74, 0.46, 0.17),
                           'Bi': (0.82, 0.71, 0.55),
                           'O': (1.00, 0.00, 0.00),
                           'S': (1.00, 1.00, 0.00),
                           'F': (0.00, 1.00, 0.00),
                           'Cl': (0.50, 1.00, 0.00),
                           'Br': (0.39, 0.15, 0.03),
                           'I': (1.00, 0.00, 1.00),
                           'Ti': (0.25, 1.75, 0.75)}
        if self.lewis is True:
            self.draw_outline = False
            self.bond_gradient = False
            self.colorbonds = False
            self.scale = 0.1
            self.ATOMS = False

    def bonds(self):
        cutoffs = self.bond_radius * covalent_radii[self.atoms.numbers]
        # ,bothways=True)
        nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
        nl.update(self.atoms)
        self.bondatoms = []
        for a in range(len(self.atoms)):
            indices, offsets = nl.get_neighbors(a)
            for a2, offset in zip(indices, offsets):
                if self.atoms[a].symbol == self.atoms[a2].symbol and self.atoms.get_distance(a, a2) <= self.double_bond:
                    bondorder = 2
                    bondoffset = self.double_bond_offset
                    if self.atoms[a].symbol == self.atoms[a2].symbol and self.atoms.get_distance(a, a2) <= self.triple_bond:
                        bondorder = 3
                        bondoffset = self.triple_bond_offset
                else:
                    bondorder = 1
                    bondoffset = [0, 0, 0]
                if (a, a2) not in self.bondatoms:
                    if self.use_bondorders is True:
                        self.bondatoms.append(
                            (a, a2, offset, bondorder, bondoffset))
                    else:
                        self.bondatoms.append((a, a2, offset, 1, [0, 0, 0]))

    def plot_bond(self, xatom1, yatom1, xatom2, yatom2, a1, a2, allx, ally, offset):
        def check_color(list1, list2):
            for n, c in enumerate(list1):
                if c != list2[n]:
                    return False
                    break
            return True
        ZORDER = min([self.atoms.positions[a1][self.view],
                     self.atoms.positions[a2][self.view]])
        colormid = (self.color1+self.color2)/2
        ca1 = self.atoms.get_chemical_symbols()[a1]
        ca2 = self.atoms.get_chemical_symbols()[a2]
        if self.draw_outline is True:
            OUTLINE = self.outline
        else:
            OUTLINE = None
        if self.bond_gradient is True:
            NUM = 61
        else:
            NUM = 3

        if np.all(offset == 0):

            if check_color(self.color1, self.color2) is True:
                self.ax.plot(allx, ally, color=self.color1, zorder=ZORDER, path_effects=[
                             self.outline], linewidth=self.bondlinewidth, solid_capstyle=self.capstyle)
            else:
                color_gradient1 = color_gradient(allx, ally, c1=self.color1, c2=self.color2, zorder=ZORDER,
                                                 linewidth=self.bondlinewidth, capstyle=self.capstyle, outline=OUTLINE, NUM=NUM)
                self.ax.add_collection(color_gradient1)
        else:
            if self.draw_outline is True:
                self.ax.plot(xatom1, yatom1, color='k', zorder=ZORDER,
                             linewidth=self.outline_bonds, solid_capstyle=self.capstyle)
                self.ax.plot(xatom2, yatom2, color='k', zorder=ZORDER,
                             linewidth=self.outline_bonds, solid_capstyle=self.capstyle)
            if ca1 == ca2:
                NUM = 3
            color_gradient1 = color_gradient(xatom1, yatom1, self.color1, colormid, zorder=ZORDER,
                                             linewidth=self.bondlinewidth, capstyle=self.capstyle, outline=None, NUM=NUM-1)
            self.ax.add_collection(color_gradient1)
            color_gradient2 = color_gradient(xatom2, yatom2, self.color2, colormid, zorder=ZORDER,
                                             linewidth=self.bondlinewidth, capstyle=self.capstyle, outline=None, NUM=NUM-1)
            self.ax.add_collection(color_gradient2)

    def plot_atoms_2D(self):
        def draw_bonds(self):
            self.bonds()
            if self.colorbonds is False:
                if self.lewis is True:
                    self.color1 = np.array((0, 0, 0))
                    self.color2 = np.array((0, 0, 0))
                else:
                    self.color1 = np.array((0.8, 0.8, 0.8))
                    self.color2 = np.array((0.8, 0.8, 0.8))

            if isinstance(self.colorbonds, str):
                self.color1 = np.array((float(self.colorbonds), float(
                    self.colorbonds), float(self.colorbonds)))
                self.color2 = self.color1
            for a1, a2, offset, bondorder, bondorderoffset in self.bondatoms:
                if np.all(offset == 0):
                    self.capstyle = 'projecting'
                    self.outline_capstyle = 'butt'
                    self.outline = mpe.withStroke(
                        linewidth=self.outline_bonds, foreground='black', capstyle='butt')
                else:
                    self.capstyle = 'projecting'
                    self.outline_capstyle = 'projecting'
                    if self.draw_outline is True:
                        self.outline = mpe.withStroke(
                            linewidth=self.outline_bonds, foreground='black', capstyle='projecting')
                    if self.draw_outline is False:
                        self.outline = None

                if self.colorbonds is True and type(self.colorbonds) is not str:
                    self.color1 = np.array(self.COLORS[a1])
                    self.color2 = np.array(self.COLORS[a2])
                if bondorder == 1:
                    atompos1 = self.new_atoms[a1].position[:-1]
                    OFFSET_REAL = np.matmul(
                        offset[:-1], np.array(self.cell_2D))
                    atompos2 = self.new_atoms[a2].position[:-1]
                    mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                    midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                    xatom1 = [mida[0], atompos1[0]]
                    yatom1 = [mida[1], atompos1[1]]
                    xatom2 = [atompos2[0], midb[0]]
                    yatom2 = [atompos2[1], midb[1]]
                    allx = [atompos1[0], atompos2[0]]
                    ally = [atompos1[1], atompos2[1]]
                    self.plot_bond(xatom1, yatom1, xatom2, yatom2,
                                   a1, a2, allx, ally, offset=offset)
                if bondorder == 2:
                    b = np.array(bondorderoffset[:-1])
                    for sign in [1, -1]:
                        atompos1 = self.new_atoms[a1].position[:-1]+b*sign
                        OFFSET_REAL = np.matmul(
                            offset[:-1], np.array(self.cell_2D))
                        atompos2 = self.new_atoms[a2].position[:-1]+b*sign
                        mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                        midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                        xatom1 = [mida[0], atompos1[0]]
                        yatom1 = [mida[1], atompos1[1]]
                        xatom2 = [atompos2[0], midb[0]]
                        yatom2 = [atompos2[1], midb[1]]
                        allx = [atompos1[0], atompos2[0]]
                        ally = [atompos1[1], atompos2[1]]
                        self.plot_bond(xatom1, yatom1, xatom2,
                                       yatom2, a1, a2, allx, ally, offset=offset)
                if bondorder == 3:
                    b = np.array(bondorderoffset[:-1])
                    for sign in [1, 0, -1]:
                        atompos1 = self.new_atoms[a1].position[:-1]+b*sign
                        OFFSET_REAL = np.matmul(
                            offset[:-1], np.array(self.cell_2D))
                        atompos2 = self.new_atoms[a2].position[:-1]+b*sign
                        mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                        midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                        xatom1 = [mida[0], atompos1[0]]
                        yatom1 = [mida[1], atompos1[1]]
                        xatom2 = [atompos2[0], midb[0]]
                        yatom2 = [atompos2[1], midb[1]]
                        allx = [atompos1[0], atompos2[0]]
                        ally = [atompos1[1], atompos2[1]]
                        self.plot_bond(xatom1, yatom1, xatom2,
                                       yatom2, a1, a2, allx, ally, offset=offset)
        self.cell_2D = [self.atoms.cell[0, :-1], self.atoms.cell[1, :-1]]
        self.real_space_shift = np.matmul(
            np.array([self.shift[0], self.shift[1]]), np.array(self.cell_2D))
        X = [0, 1, 1, 0, 0]
        Y = [0, 0, 1, 1, 0]
        self.uXs = []
        self.uYs = []
        for x, y in zip(X, Y):
            mult = np.matmul(
                np.array([x+self.shift[0], y+self.shift[1]]), np.array(self.cell_2D))
            self.uXs.append(mult[0])
            self.uYs.append(mult[1])

        if self.show_unit_cell is True:
            self.ax.plot(self.uXs, self.uYs, linestyle=self.unit_cell_linestyle, lw=self.bondlinewidth/2,
                         color='black', zorder=100, dash_capstyle='round', dash_joinstyle='round')  # zorder=10)
        self.COLORS = []
        self.sizes = []
        self.outline_bonds = self.bondlinewidth+self.bondlinewidth/2
        outline_atoms = self.bondlinewidth/5
        self.new_atoms = self.atoms.copy()
        self.new_atoms.positions = [[x+self.real_space_shift[0], y +
                                     self.real_space_shift[1], z] for x, y, z in self.atoms.positions]
        Xs = [x for x, y, z in self.new_atoms.positions]
        Ys = [y for x, y, z in self.new_atoms.positions]
        for n, a in enumerate(self.atoms.get_atomic_numbers()):
            if self.color_dict is None:
                colorval = colors.jmol_colors[a]
                if self.plot_atom_colorscaling is True:
                    colorval = colorval * \
                        self.atoms[n].position[-1] / \
                        (np.linalg.norm(np.amax(self.atoms.positions[:, 2])))
                self.COLORS.append(colorval)
            else:
                self.COLORS = [self.color_dict[s] if s in self.color_dict else colors.jmol_colors[a]
                               for s in self.atoms.get_chemical_symbols()]
            if self.lewis is True:
                self.COLORS = [(0, 0, 0)
                               for s in self.atoms.get_chemical_symbols()]
            sizeval = covalent_radii[a]*self.scale  # /np.amax(self.XARRAY)
            self.sizes.append(sizeval)
            if self.ATOMS is True:
                if self.draw_outline is True:
                    # ,zorder=10) #plot the atoms.
                    self.ax.scatter(Xs[n], Ys[n], color=self.COLORS[n], s=self.sizes[n], linewidth=outline_atoms,
                                    zorder=self.atoms.positions[n][self.view]+0.1, edgecolors='black')
                else:
                    self.ax.scatter(Xs[n], Ys[n], color=self.COLORS[n], s=self.sizes[n],
                                    zorder=self.atoms.positions[n][self.view]+0.1)
        if self.constraints is True:
            if self.atoms.constraints:
                for c in self.atoms.constraints[0].index:
                    if self.ATOMS is True:
                        self.ax.scatter(Xs[c], Ys[c], color=(0, 0, 0), marker='x', s=self.sizes[c]/(np.sqrt(
                            math.pi)), zorder=self.atoms[c].position[self.view]+0.11, alpha=0.8, linewidths=1)
        if self.show_bonds is True:
            draw_bonds(self)

    def plot_atoms_3D(self):
        def draw_bonds(self):
            self.bonds()
            # self.list_of_atoms=[a1 for a1,a2,offset,bondorder,bondorderoffset in self.bondatoms]
            # self.list_of_atoms+=[a2 for a1,a2,offset,bondorder,bondorderoffset in self.bondatoms]
            if self.draw_outline is True:
                self.outline = mpe.withStroke(
                    linewidth=self.outline_width, foreground='black', capstyle='butt')
            if self.draw_outline is False:
                self.outline = None
            if self.colorbonds is False:
                color1 = (0.8, 0.8, 0.8)
                color2 = (0.8, 0.8, 0.8)
                if self.lewis is True:
                    color1 = (0, 0, 0)
                    color2 = (0, 0, 0)
            if isinstance(self.colorbonds, str):
                color1 = (float(self.colorbonds), float(
                    self.colorbonds), float(self.colorbonds))
                color2 = color1
            for a1, a2, offset, bondorder, bondorderoffset in self.bondatoms:
                if bondorder == 1:
                    atompos1 = self.atoms[a1].position
                    OFFSET_REAL = np.matmul(offset, np.array(self.atoms.cell))
                    atompos2 = self.atoms[a2].position
                    mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                    midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                    xatom1 = [mida[0], atompos1[0]]
                    yatom1 = [mida[1], atompos1[1]]
                    zatom1 = [mida[2], atompos1[2]]
                    xatom2 = [atompos2[0], midb[0]]
                    yatom2 = [atompos2[1], midb[1]]
                    zatom2 = [atompos2[2], midb[2]]
                if bondorder == 2:
                    b = np.array(bondorderoffset)
                    for sign in [1, -1]:
                        atompos1 = self.atoms[a1].position+b*sign
                        OFFSET_REAL = np.matmul(
                            offset, np.array(self.atoms.cell))
                        atompos2 = self.atoms[a2].position+b*sign
                        mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                        midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                        zatom1 = [mida[2], atompos1[2]]
                        xatom1 = [mida[0], atompos1[0]]
                        yatom1 = [mida[1], atompos1[1]]
                        xatom2 = [atompos2[0], midb[0]]
                        yatom2 = [atompos2[1], midb[1]]
                        zatom2 = [atompos2[2], midb[2]]
                if bondorder == 3:
                    b = np.array(bondorderoffset)
                    for sign in [1, 0, -1]:
                        atompos1 = self.atoms[a1].position+b*sign
                        OFFSET_REAL = np.matmul(
                            offset, np.array(self.atoms.cell))
                        atompos2 = self.atoms[a2].position+b*sign
                        mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                        midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                        xatom1 = [mida[0], atompos1[0]]
                        yatom1 = [mida[1], atompos1[1]]
                        zatom1 = [mida[2], atompos1[2]]
                        xatom2 = [atompos2[0], midb[0]]
                        yatom2 = [atompos2[1], midb[1]]
                        zatom2 = [atompos2[2], midb[2]]

                if (type(self.colorbonds) is not float) and (self.colorbonds is True):
                    color1 = self.COLORS[a1]
                    color2 = self.COLORS[a2]
                if self.draw_outline is True:
                    self.ax.plot(xatom1, yatom1, zatom1, color=color1, zorder=atompos1[self.view], linewidth=self.bondlinewidth, path_effects=[
                                 self.outline], solid_capstyle='butt')
                    self.ax.plot(xatom2, yatom2, zatom2, color=color2, zorder=atompos2[self.view], linewidth=self.bondlinewidth, path_effects=[
                                 self.outline], solid_capstyle='butt')
                else:
                    self.ax.plot(xatom1, yatom1, zatom1, color=color1,
                                 zorder=atompos1[self.view], linewidth=self.bondlinewidth, solid_capstyle='round')
                    self.ax.plot(xatom2, yatom2, zatom2, color=color2,
                                 zorder=atompos2[self.view], linewidth=self.bondlinewidth, solid_capstyle='round')
        tags = [a.position[2] for a in self.atoms]
        self.atoms = sort(self.atoms, tags=tags)
        X = [0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
        Y = [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
        Z = [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]
        uXs = []
        uYs = []
        uZs = []
        for x, y, z in zip(X, Y, Z):
            mult = np.matmul(np.array([x, y, z]), np.array(self.atoms.cell))
            uXs.append(mult[0])
            uYs.append(mult[1])
            uZs.append(mult[2])
        if self.show_unit_cell is True:
            self.ax.plot(uXs, uYs, uZs, linestyle='dashed',
                         lw=self.bondlinewidth, color='black', zorder=0)
        Xs = [x for x, y, z in self.atoms.positions]
        Ys = [y for x, y, z in self.atoms.positions]
        Zs = [z for x, y, z in self.atoms.positions]
        self.COLORS = []
        self.sizes = []
        self.outline_bonds = self.bondlinewidth+2*self.bondlinewidth/5
        outline_atoms = self.bondlinewidth/5
        for n, a in enumerate(self.atoms.get_atomic_numbers()):
            if self.color_dict is None:
                colorval = colors.jmol_colors[a]
                if self.plot_atom_colorscaling is True:
                    colorval = colorval * \
                        self.atoms[n].position[-1] / \
                        (np.linalg.norm(np.amax(self.atoms.positions[:, 2])))
                self.COLORS.append(colorval)
            else:
                self.COLORS = [self.color_dict[s]
                               for s in self.atoms.get_chemical_symbols()]
            if self.lewis is True:
                self.COLORS = [(0, 0, 0)
                               for s in self.atoms.get_chemical_symbols()]
            sizeval = covalent_radii[a]*self.scale  # /np.amax(self.XARRAY)
            self.sizes.append(sizeval)

#            self.ax.set_xlim([uXs[0],uXs[1]])
#            self.ax.set_ylim([uYs[0],uYs[2]])
#            self.ax.set_zlim([uZs[0],uZs[6]])
        self.ax.set_box_aspect((uXs[1], uYs[2], uZs[6]))
        draw_bonds(self)
        for n, a in enumerate(self.atoms):
            if self.ATOMS is True:
                self.ax.scatter(Xs[n], Ys[n], Zs[n], color=self.COLORS[n], s=self.sizes[n],
                                zorder=a.position[self.view]+0.1, linewidth=outline_atoms, edgecolors='black', alpha=1)
        if self.ATOMS is True:
            if self.atoms.constraints and self.constraints is True:
                for c in self.atoms.constraints[0].index:
                    self.ax.scatter(Xs[c], Ys[c], Zs[c], color=(0, 0, 0), s=self.sizes[c]/(np.sqrt(
                        math.pi)), marker='x', zorder=self.atoms[c].position[self.view]+0.11, alpha=0.5)

    def plot(self):
        # if self.show is True:
        #    mpl.use('TkAgg')
        # else:
        #     mpl.use('pgf')
        if self.dimension == '2D':
            self.fig = plt.figure()
            self.ax = plt.gca()
            self.ax.set_aspect('equal')
            self.fig.patch.set_facecolor('white')
            plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
            plt.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left=False,      # ticks along the bottom edge are off
                right=False,         # ticks along the top edge are off
                labelleft=False)
            plt.axis('off')
            self.ax.set_xlim(-2, np.linalg.norm(self.atoms.cell[0])+2)
            self.ax.set_ylim(-2, np.linalg.norm(self.atoms.cell[1])+2)
            self.plot_atoms_2D()
        else:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(
                111, projection='3d', computed_zorder=False)
            self.ax.set_aspect('equal')
            self.ax.set_xlim(0, np.linalg.norm(self.atoms.cell[0]))
            self.ax.set_ylim(0, np.linalg.norm(self.atoms.cell[1]))
            self.ax.set_zlim(0, np.linalg.norm(self.atoms.cell[2]))
            self.fig.patch.set_facecolor('white')
            plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
            plt.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left=False,      # ticks along the bottom edge are off
                right=False,         # ticks along the top edge are off
                labelleft=False)
            plt.tick_params(
                axis='z',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left=False,      # ticks along the bottom edge are off
                right=False,         # ticks along the top edge are off
                labelleft=False)
            # plt.axis('off')
            self.plot_atoms_3D()
            self.ax.view_init(elev=self.elev, azim=self.azim)  # 9090
            self.ax.set_axis_off()
            self.ax.set_proj_type(self.projectiontype)
        if self.show is True:
            plt.savefig(f'{self.name}.{self.format}', transparent=True)
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(f'{self.name}.{self.format}', transparent=True)
            plt.close()
