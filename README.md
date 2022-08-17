# wmaplike

A WMAP likelihood implementation for cobaya. This likelihood is based on the original WMAP v5 likelihood implementation [(Dunkley et al., 2008)](https://arxiv.org/abs/0803.0586) as it was used in the ACT DR4 analysis [(Aiola et al., 2020)](https://arxiv.org/abs/2007.07288). This likelihood includes the high-l TT/TE/EE components, the low-l TT and beam/pointsource correction terms as they were described in the original WMAP v5 likelihood. Additionally, the high-l TB likelihood is included (although it is disabled by default). The low-l polarization likelihood is (currently) **not** included as it was unused for the original ACT DR4 analysis (see sections 6.2.1 and 6.2.2 of _Aiola et al._ for a further description).

# Installation

Clone the git with

``git clone https://github.com/HTJense/wmaplike <path/to/download>``

You can install the likelihood by running

``pip install -e <path/to/download>``

This will install the package for you. To download the WMAP DR5 data that is used for the likelihood, you need to install the likelihood with cobaya:

``cobaya-install <path/to/download> -p <cobaya/packages/path>``

This will download the WMAP data into your packages path set by your cobaya environment (this will take a while).

You can test your _WMAPLike_ installation by running `python test_wmap.py`. If all went well, the environment will give a log of the components included in WMAP, their chi-square values, and any deviations from the expected value (as reported by the original WMAP v5 likelihood with the same input data). _Note that there is currently a minor mistake in the TT beam/pointsource correction term, causing the chi-square value to deviate from the expected value by ~0.04. Because this change is so tiny, and because this component was not used for the ACT DR4 analysis either, it is not something to worry about right now, although the testing file will report on it._

# Usage

To use the likelihood, two sample .yaml file are included. The file `eval_wmap.yaml` will evaluate WMAP at the best-fitting cosmology as reported in table 2 of _Dunkley et al._ It should report on a chi-square value of 5923.52.

The file `mcmc_wmap.yaml` will run a simple MCMC chain on the likelihood. It included a prior on _tau_ (the optical depth at reionization), since there is no low-l polarization component to constrain the tau/As degeneracy otherwise.

## Settings

WMAP comes with the following settings. They are all copied, verbatim, from the original `wmap_v5` likelihood. You can check the description given in the `WMAPLike.yaml` file for descriptions.

- `use_lowl_TT` and `use_gibbs`: whether or not to use the low-ell temperature data, and whether or not to evaluate it using Dunkley et al's gibbs sampling method. Note that disabling `use_gibbs` is, as of writing, not implemented. It should refer to using the pixel likelihood, but this is not included as the gibbs sampler is faster and does not noticably affect the results.

- `use_lowl_pol`: whether or not to include the low-ell polarisation data. This is commonly replaced by a prior on tau. Should you do this, it is advised to set `temin: 24` as well, since this likelihood mimics the original WMAP's likelihood's behaviour of setting it to 2 by default.

- `use_lowl_TBEB`: whether or not to include the low-l TB and EB data. _Note that this function is currently not implemented and will raise an exception if you enable it._

- `use_highl_TT`, `use_highl_TE`, and `use_highl_TB`: whether to enable the high-l MASTER TT/TE/TB likelihoods. By default, only TT and TE are included, but TB is implemented and can be enabled if you wish.

- `use_highl_TT_beam_ptsrc`: whether or not to include the beam and pointsource correction terms in the likelihoods. By default, these dare disabled, but they are implemented, _but they do include a small mistake in the calculation._

- `use_sz`: whether or not to include the SZ template implemented in the CosmoMC implementation of the WMAP likelihood. If you do this, you can point the likelihood to a `sz_filename` (it will download the LAMBDA-provided template if it cannot find any), and you can sample over the SZ amplitude `A_sz`.

- The `ttmin`, `ttmax`, `temin` and `temax` parameters refer to the range of l-modes that should be considered by the likelihood. Note that the `temin` and `temax` parameters are also used by the B-modes.

The remainder of the settings mostly refer to names of files. You can modify these if you wish.

## The best-fitting v5 Theory

Packaged with the likelihood is the `WMAPBestFitv5` theory code that simply returns the best-fitting power spectrum that was included in the `wmap_v5` likelihood. It has one parameter, `cl_amp`, that simply scales the entire power spectrum up or down. If you want to include it, simply replace your `camb` or `classy` theory code by `WMAPBestFitv5`.
