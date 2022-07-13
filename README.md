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

You can test your _WMAPLike_ installation by running `python test_wmap.py`. If all went well, the environment will give a log of the components included in WMAP, their chi-square values, and any deviations from the expected value (as reported by the original WMAP v5 likelihood with the same input data). _Note that there is currently a minor mistake in the TT beam/pointsource correction term, causing the chi-square value to deviate from the expected value by ~0.04. Because this change is so tiny, it is not something to worry about right now, but the testing file will report on it._

# Usage

To use the likelihood, two sample .yaml file are included. The file `eval_wmap.yaml` will evaluate WMAP at the best-fitting cosmology as reported in table 2 of _Dunkley et al._ It should report on a chi-square value of 5923.52.

The file `mcmc_wmap.yaml` will run a simple MCMC chain on the likelihood. It included a prior on _tau_ (the optical depth at reionization), since there is no low-l polarization component to constrain the tau/As degeneracy otherwise.
