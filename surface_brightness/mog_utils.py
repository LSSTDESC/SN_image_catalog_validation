import numpy as np
import pandas as pd
import units
from itertools import product

def _format_extragal_catalog(galaxies, save_to_disk=None):
    # Unit conversion and column renaming
    filters = list('ugrizY')
    galaxies[['ra', 'dec']] = units.deg_to_arcsec(galaxies[['ra_true', 'dec_true']])
    galaxies['disk_to_total_ratio'] = 1.0 - galaxies['bulge_to_total_ratio_i']
    for bp in filters:
        galaxies['flux_%s' %bp] = units.mag_to_flux(galaxies['mag_true_%s_lsst' %bp].values, to_unit='nMgy')
        galaxies['flux_disk_%s' %bp] = galaxies['flux_%s' %bp].values*galaxies['disk_to_total_ratio'].values
        galaxies['flux_bulge_%s' %bp] = galaxies['flux_%s' %bp].values*galaxies['bulge_to_total_ratio_i'].values
    for component in ['disk', 'bulge']:
        galaxies['galaxy_id_%s' %component] = galaxies['galaxy_id'].values
        galaxies['e_%s' %component], galaxies['phi_%s' %component] = units.e1e2_to_ephi(e1=galaxies['ellipticity_1_%s_true' %component].values,
                   e2=galaxies['ellipticity_2_%s_true' %component].values)
        galaxies['ra_%s' %component] = galaxies['ra'].values
        galaxies['dec_%s' %component] = galaxies['dec'].values
        galaxies['size_circular_%s' %component] = (galaxies['size_minor_%s_true' %component].values*galaxies['size_%s_true' %component].values)**0.5
    # Only keep columns we'll use
    galaxies_cols_to_keep = ['galaxy_id_bulge', 'galaxy_id_disk',] #'ra', 'dec'] # 'agn', 'sprinkled', 'star']
    galaxies_cols_to_keep += ['ra_bulge', 'dec_bulge', 'size_circular_bulge', 'e_bulge', 'phi_bulge']
    galaxies_cols_to_keep += ['ra_disk', 'dec_disk', 'size_circular_disk', 'e_disk', 'phi_disk'] 
    galaxies_cols_to_keep += [prop + '_' + bp for prop, bp in product(['flux_bulge', 'flux_disk'], filters)]
    galaxies = galaxies[galaxies_cols_to_keep]
    if save_to_disk is not None:
        galaxies.to_csv(save_to_disk, index=False)
    return galaxies

def separate_bulge_disk(extragal_df):
    # Rename for convenience
    df = extragal_df
    # Separate df into bulge-related and disk-related
    bulge_df = df.filter(like='bulge', axis=1).copy()
    disk_df = df.filter(like='disk', axis=1).copy()
    # Make column schema the same across bulge and disk DataFrames (not sure if necessary)
    bulge_df.columns = [col.strip().replace('_bulge', '') for col in bulge_df.columns]
    disk_df.columns = [col.strip().replace('_disk', '') for col in disk_df.columns]
    bulge_df['sersic'] = 4
    disk_df['sersic'] = 1
    return bulge_df, disk_df, df

def sersic_to_mog(sersic_df, bulge_or_disk, flux_colname='flux_%s', filters='ugrizY'):
    from scipy.special import gammaincinv
    if bulge_or_disk=='bulge':
        # Mixture of gaussian parameters for de Vaucouleurs profile from HL13
        weights = [0.00139, 0.00941, 0.04441, 0.16162, 0.48121, 1.20357, 2.54182, 4.46441, 6.22821, 6.15393]
        stdevs = [0.00087, 0.00296, 0.00792, 0.01902, 0.04289, 0.09351, 0.20168, 0.44126, 1.01833, 2.74555]
        mog_params = {'weight': weights, 'stdev': stdevs}
        sersic_norm = gammaincinv(8, 0.5)
        gauss_norm = 40320.0*np.pi*np.exp(sersic_norm)/sersic_norm**8.0
    elif bulge_or_disk=='disk':
        # Mixture of gaussian parameters for exponential profile from HL13
        weights = [0.00077, 0.01017, 0.07313, 0.37188, 1.39727, 3.56054, 4.74340, 1.78732]
        stdevs = [0.02393, 0.06490, 0.13580, 0.25096, 0.42942, 0.69672, 1.08879, 1.67294]
        mog_params = {'weight': weights, 'stdev': stdevs}
        sersic_norm = gammaincinv(2, 0.5) # for exponential
        gauss_norm = 2.0*np.pi*np.exp(sersic_norm)/sersic_norm**2.0
    else:
        raise ValueError("Component is either bulge or disk.")
    
    mog_params_df = pd.DataFrame.from_dict(mog_params)
    # Join bulge_df and mog_params_df
    sersic_df = sersic_df.reset_index()
    sersic_df['key'] = 0
    mog_params_df['key'] = 0
    mog_df = sersic_df.merge(mog_params_df, how='left', on='key')
    mog_df = mog_df.drop('key', 1)
    mog_df['gauss_sigma'] = mog_df['size_circular']*mog_df['stdev']
    for bp in filters:
        mog_df[flux_colname %bp] = mog_df[flux_colname %bp]*mog_df['weight']/gauss_norm
    # Column sanity check
    output_cols = ['galaxy_id', 'ra', 'dec', 'e', 'phi', 'gauss_sigma',]
    output_cols += [flux_colname %bp for bp in filters]
    mog_df = mog_df[output_cols]
    return mog_df

def get_2d_gaussian(flux_array, sigma_array, ellip_array, phi_array):
    from astropy.modeling import models
    sig_sq = sigma_array**2.0
    q_sqrt = ((1.0 - ellip_array)/(1.0 + ellip_array))**0.5
    lam1 = sig_sq/q_sqrt
    lam2 = sig_sq*q_sqrt
    cos_phi = np.cos(phi_array)
    sin_phi = np.sin(phi_array)
    cov_11 = lam1*cos_phi**2.0 + lam2*sin_phi**2.0
    cov_22 = lam1*sin_phi**2.0 + lam2*cos_phi**2.0
    det = lam1*lam2
    one_norm = 1.0/(2.0*np.pi*(det)**0.5)
    assert one_norm.shape == flux_array.shape
    amp = flux_array*one_norm
    gaus_2d = models.Gaussian2D(amplitude=amp, x_mean=0.0, y_mean=0.0, x_stddev=cov_11**0.5, y_stddev=cov_22**0.5, theta=phi_array)
    return gaus_2d

def evaluate_2d_gaussian(x_array, y_array, gaus_2d_func):
    evaluated = gaus_2d_func(x_array, y_array)
    return evaluated

def sample_from_chosen_gaussian(gaus):
    mean = [0, 0]
    sig_sq = gaus['gauss_sigma']**2.0
    ellip = gaus['e']
    phi =  gaus['phi']
    q_sqrt = ((1.0 - ellip)/(1.0 + ellip))**0.5
    lam1 = sig_sq/q_sqrt
    lam2 = sig_sq*q_sqrt
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cov_11 = lam1*cos_phi**2.0 + lam2*sin_phi**2.0
    cov_22 = lam1*sin_phi**2.0 + lam2*cos_phi**2.0
    cov_12 = (lam1 - lam2)*cos_phi*sin_phi
    cov_mat = np.array([[cov_11, cov_12], [cov_12, cov_22]])
    sample = np.random.multivariate_normal(mean, cov=cov_mat)
    return sample
    