import os, sys, glob
import numpy as np
import pandas as pd
import healpy as hp
import GCRCatalogs
from GCR import GCRQuery
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
from argparse import ArgumentParser



def main(out):
    """
    Parameters
    ----------
    absolute path to output file name
    """
    # Load SN data
    sn_db = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS_DDF.db'
    engine = create_engine('sqlite:///' + sn_db)
    sn_df = pd.read_sql_table('sne_params', con=engine)
    # Note very rarely this could be a problem with 2 SN occuring in the same galaxy
    sn_df.set_index('galaxy_id', inplace=True)
    
    # Galaxies
    gals_dir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/ddf_region_galaxy_catalog/'
    fnames = glob.glob(os.path.join(gals_dir, 'gals*.hdf'))
    
    results = []
    for fname in fnames:
        gal_df = pd.read_hdf(fname)
        gal_df.set_index('galaxy_id', inplace=True)
        res = sn_df.join(gal_df).dropna()
        if len(res) > 0:
            results.append(res)
    
    joined_table = pd.concat(results)
    # write out to a csv file in the gals directory
    if out is not None:
        joined_table.to_csv(out)
    return 0

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--output_filename', help='absolute path to the output filename',
                        default=None, type=str)
    args = parser.parse_args()
    if args.output_filename == 'replace':
        gdir = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/ddf_region_galaxy_catalog/'
        out = os.path.join(gdir, 'DDF_sn_host_pairs.csv')
    else:
        out = args.output_filename
    main(out)
