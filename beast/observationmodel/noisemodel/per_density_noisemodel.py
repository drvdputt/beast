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
