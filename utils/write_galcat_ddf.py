import sys
import GCRCatalogs
from GCR import GCRQuery
import os
import argparse

## check version
#print('GCRCatalogs =', GCRCatalogs.__version__, '|' ,'GCR =', GCRCatalogs.GCR.__version__)
#for key in GCRCatalogs.available_catalogs:
#    if 'cosmo' in key:
#        print(key)

import pandas as pd


def write_hdf(healpixel, out_dir, catalog_name):
    healpix_query = GCRQuery('healpix_pixel==%d' % healpixel)
    z_filter = (lambda x: x<=1.5, 'redshift_true')
    gc = GCRCatalogs.load_catalog(catalog_name)

    x = gc.get_quantities(['galaxy_id', 'stellar_mass',
                           'stellar_mass_disk', 'stellar_mass_bulge',
                           'redshift_true', 'size_disk_true',
                           'size_minor_disk_true', 'size_minor_bulge_true',
                           'size_bulge_true', 'ra', 'dec',
                           'position_angle_true', 'ra_true', 'dec_true',
                           'morphology/diskHalfLightRadiusArcsec', 'morphology/spheroidHalfLightRadiusArcsec','morphology/positionAngle', 
                           'ellipticity_1_disk_true', 'ellipticity_2_disk_true',
                           'ellipticity_1_bulge_true', 'ellipticity_2_bulge_true',
                           'bulge_to_total_ratio_i',
                           'mag_true_u_lsst','mag_true_g_lsst', 'mag_true_r_lsst', 'mag_true_i_lsst', 'mag_true_z_lsst', 'mag_true_Y_lsst'],
                           native_filters=[healpix_query],
                           filters=[z_filter])

    df = pd.DataFrame(x)
    df.to_hdf(os.path.join(out_dir,'gals_%d_ra_dec.hdf' % healpixel),index=False, key='0')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--catalog', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    args = parser.parse_args()

    cat_config = GCRCatalogs.get_catalog_config(args.catalog)
    # healpix_list = cat_config['healpix_pixels']
    healpix_list = [8658, 8659, 8786, 8787, 8788, 8913, 8914, 8915, 8916, 9042, 9043,
       9044, 9170, 9171, 9299]
    assert os.path.isdir(args.out_dir)
    for hh in healpix_list:
        write_hdf(hh, args.out_dir, args.catalog)
