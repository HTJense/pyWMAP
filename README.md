# pyWMAP

A WMAP likelihood implementation for Cobaya. This code is a python adaptation of a subset of the publicly available WMAP v5 likelihood [(Bennet et al., 2012)](https://arxiv.org/pdf/1212.5225) retaining only data in temperature and high-l polarization. This allows to use WMAP data in combination with low-l polarization data from Planck or e.g., using a tau prior as done in the ACT DR4 analysis [(Aiola et al., 2020)](https://arxiv.org/abs/2007.07288). 
The relevant terms of the original WMAP v5 likelihood needed for this -- and that have been implemented in this adaptation of the code -- are the high-ell TT/TE/EE master components and the low-ell TT. To fully remove the TE data at low-ell the minimum multipole in TE needs to be set to ell=24 (otherwise the code uses master TE data from ell=2) and this then requires to remove the beam term in TT (the original code gives a matrix inversion error if not). The beam/pointsource correction terms and the high-l TB likelihood are also implemented in this code but disabled by default. The low-ell polarization part of the likelihood has not been implemented since better Planck data are available at those scales. 
This adaptation can be cross-checked with the original code setting to False the "use_lowl_pol" and "use_TT_beam_ptsrc" flags and setting "temin=24". This code has been coupled to Cobaya and provides identical results in likelihood/cosmological parameters compared to what can be obtained with the original code coupled to Cosmomc.

# Installation

Clone the git with

``git clone https://github.com/HTJense/wmaplike <path/to/download>``

You can install the likelihood by running

``pip install -e <path/to/download>``

This will install the package for you. To download the WMAP v5 data that is used for the likelihood, you need to install the likelihood with cobaya. From the directory `<path/to/download>`, run:

``cobaya-install wmaplike.WMAPLike -p <cobaya/packages/path>``

This will download the WMAP data into your packages path set by your cobaya environment (this will take a while).

You can test your _WMAPLike_ installation by running `python test_wmap.py`. If all went well, the environment will give a log of the components included in WMAP, their chi-square values, and any deviations from the expected value (as reported by the original WMAP v5 likelihood with the same input data). _Note that there is currently a minor mistake in the TT beam/pointsource correction term, causing the chi-square value to deviate from the expected value by ~0.04. Because this change is so tiny, and because this component should not be used, it is not something relevant, although the testing file will report on it._

# Usage

To use the likelihood, two sample .yaml file are included. The file `eval_wmap.yaml` will evaluate WMAP at the best-fitting cosmology as reported in table 2 of [Hinshaw et al. 2012](https://arxiv.org/abs/1212.5226). It should give a chi-square value of 5923.52.

The file `mcmc_wmap.yaml` will run a simple MCMC chain with the likelihood. It included a prior on _tau_, since there is no low-l polarization component to constrain the tau/As degeneracy.

## Settings

WMAP comes with the following settings. They are all copied, verbatim, from the original `wmap_v5` likelihood. You can check the description given in the `WMAPLike.yaml` file for descriptions.

- `use_lowl_TT` and `use_gibbs`: whether or not to use the low-ell temperature data, and whether or not to evaluate it using [Dunkley et al. 2008](https://arxiv.org/pdf/0803.0586)'s gibbs sampling method. Note that disabling `use_gibbs` is, as of writing, not implemented. It should refer to using the pixel likelihood, but this is not included as the gibbs sampler is faster and does not noticably affect the results.

- `use_lowl_pol`: whether or not to include the low-ell polarisation data. This is commonly replaced by a prior on tau. Should you do this, it is advised to set `temin: 24` as well, since this likelihood mimics the original WMAP's likelihood's behaviour of setting it to 2 by default.

- `use_lowl_TBEB`: whether or not to include the low-l TB and EB data. _Note that this function is currently not implemented and will raise an exception if you enable it._

- `use_highl_TT`, `use_highl_TE`, and `use_highl_TB`: whether to enable the high-l MASTER TT/TE/TB likelihoods. By default, only TT and TE are included, but TB is implemented and can be enabled if you wish.

- `use_highl_TT_beam_ptsrc`: whether or not to include the beam and pointsource correction terms in the likelihoods. By default, these dare disabled, but they are implemented, _but they do include a small mistake in the calculation._

- `use_sz`: whether or not to include the SZ template implemented in the CosmoMC implementation of the WMAP likelihood. If you do this, you can point the likelihood to a `sz_filename` (it will download the LAMBDA-provided template if it cannot find any), and you can sample over the SZ amplitude `A_sz`.

- The `ttmin`, `ttmax`, `temin` and `temax` parameters refer to the range of l-modes that should be considered by the likelihood. Note that the `temin` and `temax` parameters are also used by the B-modes.

The remainder of the settings mostly refer to names of files. You can modify these if you wish.

## The best-fitting v5 Theory

Packaged with the likelihood is the `WMAPBestFitv5` theory code that simply returns the best-fitting power spectrum that was included in the `wmap_v5` likelihood. It has one parameter, `cl_amp`, that simply scales the entire power spectrum up or down. If you want to include it, simply replace your `camb` or `classy` theory code by `WMAPBestFitv5`.
