import numpy as np
from astropy.table import Table

from ...tools.density_map import BinnedDensityMap

def split_ast_file_per_density(base_outname, astfile, binned_density_map_fname, wcs):
    """
    Splits the AST output file into parts belonging to different
    background values. Assumes that the positions are included in the
    output file, and then uses the given map to determine which bin each
    AST belongs to.

    Parameters
    ----------
    base_outname: str
        will be suffixed with the background index + '.hd5'

    astfile: str
        path to AST results

    binned_density_map: str
        path to a background map file which includes bg_bin numbers for
        each tile (the file produced by running
        pick_positions_per_background, made from a BinnedDensityMap)

    wcs: astropy WCS
        world coordinate system to convert AST positions to RA,DEC, so
        that we can find which RA,DEC tile they belong to

    """
    ast = Table.read(astfile)
    binned_density_map = BinnedDensityMap.read(binned_density_map_fname)

    xs = ast['X']
    ys = ast['Y']

    ras, decs = wcs.all_pix2world(xs, ys, 0)

    bin_foreach_ast = [binned_density_map.bin_for_position(ra, dec) for
                       ra, dec in zip(ras, decs)]

    bin_indices = np.unique(bin_foreach_ast)

    for b in bin_indices:
        idxs = np.where(bin_foreach_ast == b)
        asts_for_b = ast[idxs]
        asts_for_b.write(base_outname + '_densitybin{}.fits'.format(b))

class NoiseModelPicker:
    def __init__(self, noise_model_foreach_bin, binned_density_map_fname, wcs):
        """
        noise_model_for_each_bin should be a list of noise model
        objects, one for each density bin, in the right order.

        binned_density_map is the map that determines which bin a
        certain position belongs to, so we can pick the right one from
        noise_model_foreach_bin for each of the stars to fit.

        a wcs should be provided which will be used to convert x,y
        positions to ra,dec
        """
        #Load the map
        self.bdm = BinnedDensityMap(binned_density_map_fname)
        self.noise_model_foreach_bin = {b: noise_model_foreach_bin for b
                                        in self.bdm.bin_indices_used}
        self.wcs = wcs

    def error_and_bias_at_position(self, x, y, model_grid_index=None):
        """
        Return the noise parameters for a object based on its position
        on the image

        Parameters
        ----------
        x, y: float, float
            pixel coordinates of the object

        model_grid_index: int or list of int
            which grid indices to retrieve the noise parameters for
        """

        [ra], [dec] = self.wcs.all_pix2world(np.array([x]), np.array([y]), 0)
        b = self.bdm.bin_for_position(ra, dec)
        model = self.noise_model_foreach_bin[b]
        error = model.root.error
        bias = model.root.bias

        if model_grid_index is None:
            return error[:], bias[:]
        else:
            return error[model_grid_index], error[model_grid_index]
