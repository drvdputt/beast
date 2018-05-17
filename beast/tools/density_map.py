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

    """

    def __init__(self, bg_map, N_bg_bins=None):
        """
        bg_map: 
        Table
        of the format described above (extra columns are
        allowed)

        string
        path to file containing this table

        N_bg_bins: int
            group the tiles into background bins
        """
        if isinstance(bg_map, str):
            self.tile_data = Table.read(bg_map)
        else:
            self.tile_data = bg_map

    def tile_for_position(ra, dec):
        """
        Finds which tile a certain ra,dec fits into
        """

    def value(tile_index):
        """
        Return the background value at the given tile index
        """

    def write(fname):
        """
        Write this map to file fname (in fits format)
        """
        self.tile_data.write(fname)
        

class BinnedDensityMap(DensityMap):
    """
    Subclass which adds an extra column, which groups tiles by
    background bin. In the constructor, you can choose how many
    background bins you want.

    """
    def __init__(self, tile_data, bg_bins):
        self.tile_data = bg_map
        self.tile_data['bg_bin'] = bg_bins

    def create(self, bg_map, N_bg_bins)
        binned_density_map = DensityMap.__init__(self, bg_map)

        if N_bg_bins is None:
            bgbin_foreach_tile=np.array(range(len(self.tile_data)))

        else:
            # Create the background bins
            # [min, ., ., ., max]
            tile_bg_vals = bg['median_bg']
            min_bg = np.amin(tile_bg_vals)
            max_bg = np.amax(tile_bg_vals)
            bg_bins = np.linspace(min_bg - 0.01 * abs(min_bg),
                                  max_bg + 0.01 * abs(max_bg),
                                  N_bg_bins + 1)

            # Find which bin each tile belongs to
            # e.g. one of these numbers: 0 [1, 2, 3, 4, 5] 6
            # We have purposely chosen our bin boundaries so that no points fall
            # outside (or on the edge) of the [1,5] range
            bgbin_foreach_tile = np.digitize(tile_bg_vals, bg_bins)

        if 'bg_bin' in self.tile_data.colnames:
            print('bg_bin column already added, overwriting bg_bin column')
            self.tile_data['bg_bin']
        else:
            self.tile_data.add_column(name='bg_bin',
                                      data=bgbin_foreach_tile)

            
