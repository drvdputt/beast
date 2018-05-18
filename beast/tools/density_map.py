from astropy.table import Table, Column
import numpy as np

density_colname = 'median_bg'
bin_colname = 'bg_bin'

class DensityMap:
    """
    Class which helps with using a background map consistently, and
    allows for reliable writing and reading.

    The file format in question is the one produced by
    create_background_density_map. It should readable as an astropy
    Table, where every entry represents a tile of the background map,
    with the following tile properties as columns:

    'i_ra': right ascencion (RA) index
    'i_dec': declination (DEC) index
    'median_bg': the background value
    'min_ra': edge of tile with lowest RA
    'max_ra': edge of tile with highest RA
    'min_dec': edge of tile with lowest DEC
    'max_dec': edge of tile with highest DEC

    and the following metadata

    table.meta['ra_grid']
    table.meta['dec_grid']

    For the time being, just look at the source code of
    create_background_density_map

    """
    def __init__(self, density_map):
        """
        density_map: Table
            of the format described above (extra columns are allowed)
        --or--
        density_map: string
            path to file containing this table
        """
        if isinstance(density_map, str):
            self.tile_data = Table.read(density_map)
        else:
            self.tile_data = density_map

        self.ra_grid = self.tile_data.meta['ra_grid']
        self.dec_grid = self.tile_data.meta['dec_grid']

        self.min_i_ra = min(self.tile_data['i_ra'])
        self.max_i_ra = max(self.tile_data['i_ra'])
        self.min_i_dec = min(self.tile_data['i_dec'])
        self.max_i_dec = max(self.tile_data['i_dec'])

        # map index pairs to table rows
        self.tile_for_ij = {}
        for r in range(len(self.tile_data)):
            ij_ra_dec = (self.tile_data[r]['i_ra'], self.tile_data[r]['i_dec'])
            self.tile_for_ij[ij_ra_dec] = r

    def write(self, fname):
        """
        Write this map to file fname (in fits format)
        """
        self.tile_data.write(fname, format='hdf5', path='tile_data',
                             overwrite=True)

    def tile_for_position(self, ra, dec):
        """
        Finds which tile a certain ra,dec fits into
        """
        # Get index pair
        i_ra = np.searchsorted(self.ra_grid[:-1], side='right') - 1
        i_ra = max(i_ra, self.min_i_ra)
        i_ra = min(i_ra, self.max_i_ra)

        i_dec = np.searchsorted(self.dec_grid[:-1], side='right') - 1
        i_dec = max(i_dec, self.min_i_ra)
        i_dec = min(i_dec, self.max_i_ra)

        # Use index pair to row index map
        return self.tile_for_ij[(i_ra, i_dec)]

    def ra_dec_mins(self):
        return self.tile_data['ra_min'], self.tile_data['dec_min']

    def ra_dec_deltas(self):
        return (self.tile_data['ra_max'] - self.tile_data['ra_min'],
                self.tile_data['dec_max'] - self.tile_data['dec_min'])

    def value(self, tile_index):
        """
        Return the background value at the given tile index
        """
        return self.tile_data[density_colname][tile_index]


class BinnedDensityMap(DensityMap):
    """
    Subclass which adds an extra column, which groups tiles by
    background bin. It is recommended to not use the constructor
    directly. In the 'create' function, you can choose how many
    background bins you want, while the 'read' function reads an
    existing binned density map from file.

    """

    def __init__(self, tile_data, bins):
        self.tile_data = tile_data
        self.tile_data[bin_colname] = bins
        self.bin_indices_used = np.sort(np.unique(bins))

    def create(density_map, N_bins=None):
        """
        Creates a binned density map from a DensityMap file, or from an
        astropy table loaded from it. The tiles are grouped into
        N_bins density bins.

        If N_bins is none, each tile is treated as a separate bin.
        """
        binned_density_map = DensityMap(density_map)

        if N_bins is None:
            bin_foreach_tile = np.array(
                range(len(binned_density_map.tile_data)))

        else:
            # Create the background bins
            # [min, ., ., ., max]
            tile_densities = binned_density_map.tile_data[density_colname]
            min_density = np.amin(tile_densities)
            max_density = np.amax(tile_densities)
            bins = np.linspace(min_density - 0.01 * abs(min_density),
                               max_density + 0.01 * abs(max_density),
                               N_bins + 1)

            # Find which bin each tile belongs to
            # e.g. one of these numbers: 0 [1, 2, 3, 4, 5] 6
            # We have purposely chosen our bin boundaries so that no points fall
            # outside (or on the edge) of the [1,5] range
            bin_foreach_tile = np.digitize(tile_densities, bins)

        if bin_colname in binned_density_map.tile_data.colnames:
            print('{} column already there, overwriting it.'.format(bin_colname))
            binned_density_map.tile_data[bin_colname] = bin_foreach_tile
        else:
            c = Column(name=bin_colname, data=bin_foreach_tile)
            binned_density_map.tile_data.add_column(c)

    def read(density_map_fname):
        binned_density_map = DensityMap(density_map_fname)
        if bin_colname not in binned_density_map.tile_data.colnames:
            raise Exception('{} column not yet calculated. '
                            'Please use \'read\' function instead'.format(bin_colname))

    def bin_for_position(self, ra, dec):
        """
        Finds which background bin a certain ra,dec fits into, and
        returns its index.
        """
        t = self.tile_for_position(ra, dec)
        return self.tile_data[bin_colname][t]

    def value_foreach_tile(self):
        return self.tile_data[density_colname]

    def bin_foreach_tile(self):
        return self.tile_data[bin_colname]

    def tiles_foreach_bin(self):
        """
        Invert the above function. Result is a list of lists, each
        nested list representing a set of tiles that belong to the same
        bin.
        """
        return [ # Find where bgbin_foreach_tile is equal to b
            np.nonzero(self.bgbin_foreach_tile == b)[0]
            # for all the active bin indices
            for b in self.bin_indices_used]
