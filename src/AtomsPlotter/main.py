from ase.io import read
import numpy as np
from ase.build import make_supercell
from .atoms_plotter import atoms_plotter


def main():
    import argparse

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n'):
            return False
        else:
            # raise argparse.ArgumentTypeError('Boolean value expected.')
            return v
    parser = argparse.ArgumentParser(
        description='Plotter using utility of ASE, in matplotlib')
    parser.add_argument(
        'atoms', type=str, help='atoms object, readable by ASE .poscar, .vasp, .xyz, ...')
    parser.add_argument('-r', '--repeat', default=[
                        1, 1, 1], help='repetition of the unit cell in [x,y] direction', type=int, nargs=3)
    parser.add_argument('-b', '--bondlinewidth',
                        help='line width of the bonds', default=1, type=float)
    parser.add_argument(
        '-B', '--bonds', help='plot the atoms but not the bonds', default=True, type=str2bool)
    # parser.add_argument('-u','--unitcell',help='decide if you want to plot the atoms but not the unit cell',type=str2bool,default=True)
    parser.add_argument(
        '-u', '--unitcell', help='plot unit cell. Buggy atoms are rotated', action='store_true')
    parser.add_argument(
        '-S', '--show', help='decide if you only want to save or show the figure', action='store_true')
    parser.add_argument('-N', '--name', default='filename',
                        help='name for the saved plot if show == False', type=str)
    parser.add_argument('-s', '--scale', default=100,
                        help='empiric scaling value for the atoms in the plot', type=int)
    parser.add_argument(
        '-c', '--colorbonds', help='use the colors of the atoms for the bonds; or give value between 1 and 0; or False', type=str2bool, default=True)
    parser.add_argument('-O', '--outline',
                        help='draw outline', action='store_true')
    parser.add_argument(
        '-j', '--jmol', help='use the jmol colors instead of default-dict', action='store_true')
    parser.add_argument('-D', '--dimension', help='plotting in 2 or 3D',
                        choices=['2D', '3D'], default='2D')
    parser.add_argument(
        '--elev', help='elevate angle; currently zorder is only correct for 90', default=90, type=float)
    parser.add_argument(
        '--azim', help='azimuth angle; currently zorder is only correct for 90', default=90, type=float)
    parser.add_argument('--unit_cell_linestyle', default='dashed', type=str)
    parser.add_argument('-C', '--cutoff_radius', default=1,
                        help='cutoff radius for drawing the bonds', type=float)
    parser.add_argument(
        '--no_bondorders', help='call if you don\'t want double or triple bonds', action='store_false')
    parser.add_argument('--triple_bonds', default=1.34,
                        help='distance below which bonds are read as triple bonds', type=float)
    parser.add_argument('--double_bonds', default=1.4,
                        help='distance below which bonds are read as double bonds', type=float)
    parser.add_argument('--double_bond_offset', default=0.1,
                        help='offset b will be used for plotting double bonds [b,b,0]', type=float)
    parser.add_argument('--triple_bond_offset', default=0.1,
                        help='offset t will be used for plotting triple bonds [t,t,0]', type=float)
    parser.add_argument(
        '--bond_gradient', help='use gradient for coloring the bonds', type=str2bool, default=True)
    parser.add_argument(
        '--constraints', help='draw constraints', type=str2bool, default=True)
    parser.add_argument(
        '--rot', help='rotation of the atoms in [x,y,z] direction', type=float, nargs=3)
    parser.add_argument(
        '-v', '--view', help="view direction for z-order, use 0,1,2 for x,y,z", default=2, type=int)
    parser.add_argument(
        '-l', '--plot_lewis', help='plot lewis structure instead of colored', action='store_true')
    a = parser.parse_args()
    atoms = read(a.atoms)
    plotter = atoms_plotter(atoms=atoms,
                            show=a.show,
                            lewis=a.plot_lewis,
                            show_unit_cell=a.unitcell,
                            repeat=a.repeat,
                            bondlinewidth=a.bondlinewidth,
                            bond_cutoff=a.cutoff_radius,
                            colorbonds=a.colorbonds,
                            draw_outline=a.outline,
                            constraints=a.constraints,
                            scale=a.scale,
                            name=a.name,
                            triple_bond=a.triple_bonds,
                            double_bond=a.double_bonds,
                            double_bond_offset=a.double_bond_offset,
                            triple_bond_offset=a.triple_bond_offset,
                            bond_gradient=a.bond_gradient,
                            dimension=a.dimension,
                            view=a.view,
                            )
    if a.rot:
        plotter.atoms.wrap(pretty_translation=True)
        plotter.atoms.rotate(a.rot[0], 'x', center='COM', rotate_cell=True)
        plotter.atoms.rotate(a.rot[1], 'y', center='COM', rotate_cell=True)
        plotter.atoms.rotate(a.rot[2], 'z', center='COM', rotate_cell=True)
    if np.any(a.repeat) != 1:
        plotter.atoms = make_supercell(plotter.atoms, [[int(a.repeat[0]), 0, 0], [
                                       0, int(a.repeat[1]), 0], [0, 0, int(a.repeat[2])]])
    if a.name == 'filename':
        plotter.name = a.atoms.split(".")[0]
    else:
        plotter.name = a.name

    scaling = (np.linalg.norm(
        plotter.atoms.cell[0])+np.linalg.norm(plotter.atoms.cell[2]))/3/20
    plotter.scale = a.scale/scaling/float(a.repeat[0])  # 3D 50
    plotter.bondlinewidth = 1.5*a.bondlinewidth/scaling/float(a.repeat[0])
    if a.no_bondorders is False:
        plotter.use_bondorders = False
    else:
        plotter.double_bonds = a.double_bonds
        plotter.double_bond_offset = [
            a.double_bond_offset, a.double_bond_offset, 0]
        plotter.triple_bonds = a.triple_bonds
        plotter.triple_bond_offset = [
            a.triple_bond_offset, a.triple_bond_offset, 0]
    if a.jmol is True:
        plotter.color_dict = None
    # projection

    plotter.azim = a.azim
    plotter.elev = a.elev
    plotter.plot()


if __name__ == "__main__":
    main()
