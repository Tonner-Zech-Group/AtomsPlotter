import math

import ase
import matplotlib.patheffects as mpe
import matplotlib.pyplot as plt
import numpy as np
from ase.build import sort
from ase.data import colors, covalent_radii
from ase.neighborlist import NeighborList

from .utils import color_gradient

# Stick width of a default bond (bondlinewidth=1) in Angstrom, used when
# sizes are derived from the structure (scale=None).
BOND_WIDTH_ANG = 0.15
# Extra margin in Angstrom around non-periodic structures.
FRAME_MARGIN_ANG = 0.25


class atoms_plotter():
    def __init__(self,
                 atoms=None,
                 ATOMS=True,
                 view=2,
                 repeat=[1, 1],
                 show_unit_cell=False,
                 bond_cutoff=1,
                 unit_cell_linestyle='dashed',
                 # None: derive marker sizes from the covalent radii in
                 # real-space units (uniform across structures). A number
                 # restores the legacy behaviour s = covalent_radius * scale
                 # in points^2.
                 scale=None,
                 atom_scale=0.5,
                 inches_per_angstrom=0.4,
                 min_figsize=3.0,
                 max_figsize=16.0,
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
                 plot_atom_colorscaling=False,
                 # None: orient non-periodic structures along their
                 # principal axes (largest extent on x, then y) so the
                 # default view shows the molecule's plane instead of an
                 # arbitrary projection. True/False forces it on/off.
                 auto_orient=None,
                 format='svg'):

        self.format = format
        self.atoms = atoms if atoms is not None else ase.Atoms()
        self.ATOMS = ATOMS

        self.view = view
        self.lewis = lewis
        self.repeat = repeat
        self.plot_atom_colorscaling = plot_atom_colorscaling
        self.show_unit_cell = show_unit_cell
        self.bond_cutoff = bond_cutoff
        self.unit_cell_linestyle = unit_cell_linestyle
        self.scale = scale
        self.atom_scale = atom_scale
        self.inches_per_angstrom = inches_per_angstrom
        self.min_figsize = min_figsize
        self.max_figsize = max_figsize
        self.show = show
        self.show_bonds = show_bonds
        self.use_bondorders = use_bondorders
        self.bondatoms = []
        self.double_bond_offset = double_bond_offset
        self.triple_bond_offset = triple_bond_offset
        self.outline_width = outline_width
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
        self.auto_orient = auto_orient
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
                           'Ti': (0.25, 1.75, 0.75),
                           'Au': (1.00, 0.84, 0.00)}
        if self.lewis is True:
            self.draw_outline = False
            self.bond_gradient = False
            self.colorbonds = False
            self.ATOMS = False

    def _orient_principal_axes(self):
        """Rotate the atoms so their principal axes align with x, y, z.

        The largest spatial extent ends up on x and the smallest on z,
        which makes the default x-y projection show e.g. a planar
        molecule face-on and a linear molecule lengthwise.
        """
        positions = self.atoms.positions
        if len(self.atoms) < 2:
            return
        centered = positions - positions.mean(axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(centered.T @ centered)
        # eigh returns ascending eigenvalues; we want largest axis first
        eigenvectors = eigenvectors[:, ::-1]
        if np.linalg.det(eigenvectors) < 0:
            eigenvectors[:, 2] *= -1
        self.atoms.positions = centered @ eigenvectors

    def _drawn_radii(self):
        """Radii (Angstrom) the atoms are drawn with, used for framing."""
        fraction = self.atom_scale if self.scale is None else 0.5
        return covalent_radii[self.atoms.numbers] * fraction

    def check_pbc(self):
        """Compute the display frame; returns True if fully periodic.

        Sets self.frame = [[xmin, xmax], [ymin, ymax], [zmin, zmax]].
        Fully periodic: the frame spans the cell, extended if atoms stick
        out of it. Otherwise: the frame is the bounding box of the atoms
        padded by their drawn radii, i.e. the structure is centered in
        the figure no matter where it sits in space. Non-periodic
        structures additionally get an axis-aligned bounding-box cell so
        downstream cell math (bond offsets, 3D aspect) has something to
        work with.
        """
        pbc = self.atoms.pbc
        positions = self.atoms.positions
        radii = self._drawn_radii() + FRAME_MARGIN_ANG
        lo = (positions - radii[:, None]).min(axis=0)
        hi = (positions + radii[:, None]).max(axis=0)
        if pbc.any():
            corners = np.array([[i, j, k]
                                for i in (0, 1)
                                for j in (0, 1)
                                for k in (0, 1)]) @ np.asarray(self.atoms.cell)
            lo = np.minimum(lo, corners.min(axis=0))
            hi = np.maximum(hi, corners.max(axis=0))
        else:
            self.atoms.cell = np.diag(hi - lo)
        self.frame = [[lo[i], hi[i]] for i in range(3)]
        self.frame_x, self.frame_y, self.frame_z = self.frame
        return bool(pbc.all())

    def bonds(self):
        cutoffs = self.bond_cutoff * covalent_radii[self.atoms.numbers]
        nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
        nl.update(self.atoms)
        cell = np.asarray(self.atoms.cell)
        self.bondatoms = []
        for a in range(len(self.atoms)):
            indices, offsets = nl.get_neighbors(a)
            for a2, offset in zip(indices, offsets):
                distance = np.linalg.norm(
                    self.atoms.positions[a2] + offset @ cell
                    - self.atoms.positions[a])
                if (self.atoms[a].symbol == self.atoms[a2].symbol
                        and distance <= self.double_bond):
                    if distance <= self.triple_bond:
                        bondorder = 3
                        bondoffset = self.triple_bond_offset
                    else:
                        bondorder = 2
                        bondoffset = self.double_bond_offset
                else:
                    bondorder = 1
                    bondoffset = [0, 0, 0]
                if self.use_bondorders is True:
                    self.bondatoms.append(
                        (a, a2, offset, bondorder, bondoffset))
                else:
                    self.bondatoms.append((a, a2, offset, 1, [0, 0, 0]))

    def _atom_colors(self):
        """Per-atom colors honoring lewis mode, color_dict and jmol fallback."""
        numbers = self.atoms.get_atomic_numbers()
        symbols = self.atoms.get_chemical_symbols()
        if self.lewis is True:
            return [(0, 0, 0) for _ in symbols]
        if self.color_dict is None:
            COLORS = [colors.jmol_colors[z] for z in numbers]
            if self.plot_atom_colorscaling is True:
                # Scale each color by the atom's z fraction. Skip flat
                # structures (z extent 0) to avoid dividing by zero.
                z = self.atoms.positions[:, 2]
                z_norm = np.abs(z).max()
                if z_norm > 0:
                    COLORS = [np.clip(np.asarray(c) * z[i] / z_norm, 0, 1)
                              for i, c in enumerate(COLORS)]
            return COLORS
        return [self.color_dict[s] if s in self.color_dict
                else colors.jmol_colors[numbers[i]]
                for i, s in enumerate(symbols)]

    def _marker_sizes(self):
        """Scatter sizes (points^2). Real-space based unless legacy scale."""
        numbers = self.atoms.get_atomic_numbers()
        if self.scale is not None:
            return [covalent_radii[z] * self.scale for z in numbers]
        return [(2 * covalent_radii[z] * self.atom_scale
                 * self._pt_per_ang) ** 2 for z in numbers]

    def plot_bond(self, xatom1, yatom1, xatom2, yatom2, a1, a2, allx, ally, offset):
        ZORDER = min([self.atoms.positions[a1][self.view],
                     self.atoms.positions[a2][self.view]])
        colormid = (self.color1 + self.color2) / 2
        ca1 = self.atoms.get_chemical_symbols()[a1]
        ca2 = self.atoms.get_chemical_symbols()[a2]
        OUTLINE = self.outline if self.draw_outline is True else None
        NUM = 61 if self.bond_gradient is True else 3

        if np.all(offset == 0):
            if np.allclose(self.color1, self.color2):
                path_effects = [OUTLINE] if OUTLINE is not None else None
                self.ax.plot(allx, ally, color=self.color1, zorder=ZORDER,
                             path_effects=path_effects,
                             linewidth=self._bond_lw,
                             solid_capstyle=self.capstyle)
            else:
                gradient = color_gradient(
                    allx, ally, c1=self.color1, c2=self.color2, zorder=ZORDER,
                    linewidth=self._bond_lw, capstyle=self.capstyle,
                    outline=OUTLINE, NUM=NUM)
                self.ax.add_collection(gradient)
        else:
            if self.draw_outline is True:
                self.ax.plot(xatom1, yatom1, color='k', zorder=ZORDER,
                             linewidth=self.outline_bonds,
                             solid_capstyle=self.capstyle)
                self.ax.plot(xatom2, yatom2, color='k', zorder=ZORDER,
                             linewidth=self.outline_bonds,
                             solid_capstyle=self.capstyle)
            if ca1 == ca2:
                NUM = 3
            gradient1 = color_gradient(
                xatom1, yatom1, self.color1, colormid, zorder=ZORDER,
                linewidth=self._bond_lw, capstyle=self.capstyle,
                outline=None, NUM=NUM - 1)
            self.ax.add_collection(gradient1)
            gradient2 = color_gradient(
                xatom2, yatom2, self.color2, colormid, zorder=ZORDER,
                linewidth=self._bond_lw, capstyle=self.capstyle,
                outline=None, NUM=NUM - 1)
            self.ax.add_collection(gradient2)

    def plot_atoms_2D(self):
        def draw_bonds(self):
            self.bonds()
            if self.colorbonds is False:
                if self.lewis is True:
                    self.color1 = np.array((0, 0, 0))
                else:
                    self.color1 = np.array((0.8, 0.8, 0.8))
                self.color2 = self.color1
            if isinstance(self.colorbonds, str):
                self.color1 = np.array([float(self.colorbonds)] * 3)
                self.color2 = self.color1
            for a1, a2, offset, bondorder, bondorderoffset in self.bondatoms:
                self.capstyle = 'projecting'
                if np.all(offset == 0):
                    capstyle = 'butt'
                else:
                    capstyle = 'projecting'
                if self.draw_outline is True:
                    self.outline = mpe.withStroke(
                        linewidth=self.outline_bonds, foreground='black',
                        capstyle=capstyle)
                else:
                    self.outline = None

                if self.colorbonds is True:
                    self.color1 = np.array(self.COLORS[a1])
                    self.color2 = np.array(self.COLORS[a2])
                b = np.array(bondorderoffset[:-1])
                signs = {1: [0], 2: [1, -1], 3: [1, 0, -1]}[bondorder]
                for sign in signs:
                    atompos1 = self.new_atoms[a1].position[:-1] + b * sign
                    atompos2 = self.new_atoms[a2].position[:-1] + b * sign
                    OFFSET_REAL = np.matmul(
                        offset[:-1], np.array(self.cell_2D))
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
        self.cell_2D = [self.atoms.cell[0, :-1], self.atoms.cell[1, :-1]]
        self.real_space_shift = np.matmul(
            np.array([self.shift[0], self.shift[1]]), np.array(self.cell_2D))
        X = [0, 1, 1, 0, 0]
        Y = [0, 0, 1, 1, 0]
        self.uXs = []
        self.uYs = []
        for x, y in zip(X, Y):
            mult = np.matmul(
                np.array([x + self.shift[0], y + self.shift[1]]),
                np.array(self.cell_2D))
            self.uXs.append(mult[0])
            self.uYs.append(mult[1])

        if self.show_unit_cell is True:
            self.ax.plot(self.uXs, self.uYs,
                         linestyle=self.unit_cell_linestyle,
                         lw=self._bond_lw / 2, color='black', zorder=100,
                         dash_capstyle='round', dash_joinstyle='round')
        self.outline_bonds = self._bond_lw * 1.5
        outline_atoms = self._bond_lw / 5
        self.new_atoms = self.atoms.copy()
        self.new_atoms.positions = self.atoms.positions + np.append(
            self.real_space_shift, 0)
        Xs = self.new_atoms.positions[:, 0]
        Ys = self.new_atoms.positions[:, 1]
        self.COLORS = self._atom_colors()
        self.sizes = self._marker_sizes()
        if self.ATOMS is True:
            for n in range(len(self.atoms)):
                if self.draw_outline is True:
                    self.ax.scatter(Xs[n], Ys[n], color=self.COLORS[n],
                                    s=self.sizes[n], linewidth=outline_atoms,
                                    zorder=self.atoms.positions[n][self.view] + 0.1,
                                    edgecolors='black')
                else:
                    self.ax.scatter(Xs[n], Ys[n], color=self.COLORS[n],
                                    s=self.sizes[n],
                                    zorder=self.atoms.positions[n][self.view] + 0.1)
            if self.constraints is True and self.atoms.constraints:
                for c in self.atoms.constraints[0].index:
                    self.ax.scatter(Xs[c], Ys[c], color=(0, 0, 0), marker='x',
                                    s=self.sizes[c] / np.sqrt(math.pi),
                                    zorder=self.atoms[c].position[self.view] + 0.11,
                                    alpha=0.8, linewidths=1)
        if self.show_bonds is True:
            draw_bonds(self)

    def _view_direction(self):
        """Unit vector pointing from the structure towards the camera."""
        elev = math.radians(self.elev)
        azim = math.radians(self.azim)
        return np.array([math.cos(elev) * math.cos(azim),
                         math.cos(elev) * math.sin(azim),
                         math.sin(elev)])

    def plot_atoms_3D(self):
        # depth along the actual camera direction, so occlusion is
        # correct for any elev/azim (not just the top-down view)
        view_dir = self._view_direction()

        def depth(position):
            return float(np.dot(position, view_dir))

        def draw_bonds(self):
            self.bonds()
            if self.draw_outline is True:
                self.outline = mpe.withStroke(
                    linewidth=max(self.outline_width, self._bond_lw * 1.5),
                    foreground='black', capstyle='butt')
            else:
                self.outline = None
            color1 = color2 = None
            if self.colorbonds is False:
                color1 = (0, 0, 0) if self.lewis is True else (0.8, 0.8, 0.8)
                color2 = color1
            if isinstance(self.colorbonds, str):
                color1 = (float(self.colorbonds),) * 3
                color2 = color1
            for a1, a2, offset, bondorder, bondorderoffset in self.bondatoms:
                if self.colorbonds is True:
                    color1 = self.COLORS[a1]
                    color2 = self.COLORS[a2]
                b = np.array(bondorderoffset)
                signs = {1: [0], 2: [1, -1], 3: [1, 0, -1]}[bondorder]
                for sign in signs:
                    atompos1 = self.atoms[a1].position + b * sign
                    atompos2 = self.atoms[a2].position + b * sign
                    OFFSET_REAL = np.matmul(offset, np.array(self.atoms.cell))
                    mida = 0.5 * (atompos1 + atompos2 + OFFSET_REAL)
                    midb = 0.5 * (atompos1 + atompos2 - OFFSET_REAL)
                    xatom1 = [mida[0], atompos1[0]]
                    yatom1 = [mida[1], atompos1[1]]
                    zatom1 = [mida[2], atompos1[2]]
                    xatom2 = [atompos2[0], midb[0]]
                    yatom2 = [atompos2[1], midb[1]]
                    zatom2 = [atompos2[2], midb[2]]
                    if self.draw_outline is True:
                        self.ax.plot(xatom1, yatom1, zatom1, color=color1,
                                     zorder=depth(atompos1),
                                     linewidth=self._bond_lw,
                                     path_effects=[self.outline],
                                     solid_capstyle='butt')
                        self.ax.plot(xatom2, yatom2, zatom2, color=color2,
                                     zorder=depth(atompos2),
                                     linewidth=self._bond_lw,
                                     path_effects=[self.outline],
                                     solid_capstyle='butt')
                    else:
                        self.ax.plot(xatom1, yatom1, zatom1, color=color1,
                                     zorder=depth(atompos1),
                                     linewidth=self._bond_lw,
                                     solid_capstyle='round')
                        self.ax.plot(xatom2, yatom2, zatom2, color=color2,
                                     zorder=depth(atompos2),
                                     linewidth=self._bond_lw,
                                     solid_capstyle='round')
        tags = [depth(a.position) for a in self.atoms]
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
            self.ax.plot(uXs, uYs, uZs, linestyle=self.unit_cell_linestyle,
                         lw=self._bond_lw, color='black', zorder=0)
        Xs = self.atoms.positions[:, 0]
        Ys = self.atoms.positions[:, 1]
        Zs = self.atoms.positions[:, 2]
        self.outline_bonds = self._bond_lw * 1.5
        outline_atoms = self._bond_lw / 5
        self.COLORS = self._atom_colors()
        self.sizes = self._marker_sizes()
        draw_bonds(self)
        if self.ATOMS is True:
            for n, a in enumerate(self.atoms):
                self.ax.scatter(Xs[n], Ys[n], Zs[n], color=self.COLORS[n],
                                s=self.sizes[n],
                                zorder=depth(a.position) + 0.1,
                                linewidth=outline_atoms, edgecolors='black',
                                alpha=1)
            if self.atoms.constraints and self.constraints is True:
                for c in self.atoms.constraints[0].index:
                    self.ax.scatter(Xs[c], Ys[c], Zs[c], color=(0, 0, 0),
                                    s=self.sizes[c] / np.sqrt(math.pi),
                                    marker='x',
                                    zorder=depth(self.atoms[c].position) + 0.11,
                                    alpha=0.5)

    def _figsize(self, width_ang, height_ang):
        """Figure size in inches for a frame of the given real-space size.

        The size is proportional to the structure (inches_per_angstrom),
        clamped so the longer edge stays within [min_figsize, max_figsize].
        """
        width = width_ang * self.inches_per_angstrom
        height = height_ang * self.inches_per_angstrom
        longest = max(width, height)
        if longest > self.max_figsize:
            factor = self.max_figsize / longest
        elif longest < self.min_figsize:
            factor = self.min_figsize / longest
        else:
            factor = 1.0
        return width * factor, height * factor

    def plot(self):
        original_atoms = self.atoms
        self.atoms = self.atoms.copy()
        try:
            repeat = list(self.repeat) + [1] * (3 - len(self.repeat))
            if self.atoms.pbc.any() and any(r != 1 for r in repeat):
                self.atoms = self.atoms.repeat(
                    tuple(int(r) for r in repeat[:3]))
            orient = self.auto_orient
            if orient is None:
                orient = not self.atoms.pbc.any()
            if orient:
                self._orient_principal_axes()
            self.check_pbc()
            self.frame_x, self.frame_y, self.frame_z = self.frame
            spans = [hi - lo for lo, hi in self.frame]

            if self.dimension == '2D':
                fig_w, fig_h = self._figsize(spans[0], spans[1])
                # Points available per Angstrom: the axes fill the whole
                # figure, so atom/bond sizes in points follow real space.
                self._pt_per_ang = 72 * fig_w / spans[0]
                if self.scale is None:
                    self._bond_lw = (self.bondlinewidth * BOND_WIDTH_ANG
                                     * self._pt_per_ang)
                else:
                    self._bond_lw = self.bondlinewidth
                self.fig = plt.figure(figsize=(fig_w, fig_h))
                self.ax = self.fig.add_axes([0, 0, 1, 1])
                self.ax.set_aspect('equal')
                self.ax.set_axis_off()
                self.fig.patch.set_facecolor('white')
                self.ax.set_xlim(self.frame[0][0], self.frame[0][1])
                self.ax.set_ylim(self.frame[1][0], self.frame[1][1])
                self.plot_atoms_2D()
            else:
                size = max(self._figsize(max(spans), max(spans)))
                self._pt_per_ang = 72 * size / max(spans)
                if self.scale is None:
                    self._bond_lw = (self.bondlinewidth * BOND_WIDTH_ANG
                                     * self._pt_per_ang)
                else:
                    self._bond_lw = self.bondlinewidth
                self.fig = plt.figure(figsize=(size, size))
                self.ax = self.fig.add_subplot(
                    111, projection='3d', computed_zorder=False)
                self.fig.patch.set_facecolor('white')
                self.ax.set_xlim(self.frame[0][0], self.frame[0][1])
                self.ax.set_ylim(self.frame[1][0], self.frame[1][1])
                self.ax.set_zlim(self.frame[2][0], self.frame[2][1])
                # zoom > 1 trims the large default margins of 3D axes
                self.ax.set_box_aspect(spans, zoom=1.3)
                self.plot_atoms_3D()
                self.ax.view_init(elev=self.elev, azim=self.azim)
                self.ax.set_axis_off()
                self.ax.set_proj_type(self.projectiontype)
            if self.show is True:
                plt.show()
            else:
                plt.savefig(f'{self.name}.{self.format}', transparent=True)
                plt.close()
        finally:
            self.atoms = original_atoms
