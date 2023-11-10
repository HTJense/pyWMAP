"""
.. module:: wmaplike

:Synopsis: A python implementation of the WMAP likelihood. See Dunkley et al.
           (2008) for a full description.
:Authors: Hidde Jense, Umberto Natale.
"""
import os
import numpy as np

from typing import Sequence
from scipy import linalg

from astropy.io import fits
from cobaya import Theory
from cobaya.likelihoods.base_classes import InstallableLikelihood
from cobaya.log import LoggedError


class WMAPLike(InstallableLikelihood):
    install_options = {
        "download_url": "https://lambda.gsfc.nasa.gov/data/map/dr5/dcp/\
                         likelihood/wmap_likelihood_full_v5.tar.gz",
        "data_path": "WMAP/dr5"
    }

    def initialize(self) -> None:
        data_file_path = os.path.normpath(getattr(self, "path", None) or
                                          os.path.join(self.packages_path,
                                                       "data"))
        self.data_folder = os.path.join(data_file_path, self.data_folder)

        if not os.path.exists(self.data_folder):
            raise LoggedError(self.log, f"No data folder found at \
                                          [{self.data_folder}].")

        self.sz_file = os.path.join(self.data_folder, self.sz_filename)

        if self.use_sz:
            self.setup_for_sz()

        self.tt_lmax = 1200
        self.te_lmax = 800

        self.TEEEBB_pixlike_logdet_offset = 16078.083180
        self.TE_logdet_offset = 3584.277805
        self.TB_logdet_offset = 3598.152208

        self.use_cl = set()

        if self.use_lowl_TT:
            """ Initialize the low l TT likelihood. """
            self.use_cl.update("tt")

            if self.use_gibbs:
                self.setup_for_tt_gibbs()
            else:
                self.setup_for_tt_exact()

            self.tt_highl_start = self.lowl_max + 1
        else:
            self.tt_highl_start = self.ttmin

        if self.use_lowl_pol:
            self.use_cl.update("te")
            self.use_cl.update("ee")
            self.use_cl.update("bb")

            self.te_highl_start = 24

            if self.use_lowl_TBEB:
                self.use_cl.update("tb")
                self.use_cl.update("eb")

                self.setup_for_lowl_te_tb_ee_bb_eb()
            else:
                self.setup_for_lowl_te_ee_bb()
        else:
            self.te_highl_start = self.temin

            if self.te_highl_start < 24:
                self.log.warn("Note: you have disabled the low-ell \
                               polarisation functions, but you did not \
                               include the data. If you intend to replace the \
                               low-ell polarisation data by a tau prior, you \
                               should include `temin: 24` as a setting.")

        if self.use_highl_TT or self.use_highl_TE:
            self.use_cl.update("tt")
            self.use_cl.update("te")

            self.setup_for_tt_te_covariance()

        if self.use_highl_TT_beam_ptsrc:
            self.use_cl.update("tt")

            self.setup_for_tt_beam_ptsrc(2, self.tt_lmax)

        if self.use_highl_TB:
            self.use_cl.update("tt")
            self.use_cl.update("tb")
            self.use_cl.update("bb")

            self.setup_for_tt_te_covariance()

    def setup_for_sz(self) -> None:
        if not os.path.exists(self.sz_file):
            self.log.info(f"Could not find the SZ spectrum file. Will try to \
                            download it to [{self.sz_file}].")

            if not os.path.exists(os.path.dirname(self.sz_file)):
                os.makedirs(os.path.dirname(self.sz_file))
                self.log.debug(f"Created new directory \
                                 {os.path.dirname(self.sz_file)}.")

            import requests

            url = "https://lambda.gsfc.nasa.gov/data/map/dr5/dcp/sz_spectra/\
                   wmap_sz_spectrum_61GHz_v5.txt"

            req = requests.get(url, stream=True)
            if req.ok:
                self.log.info("Downloading...")

                with open(self.sz_file, "wb") as fp:
                    for chunk in req.iter_content(chunk_size=1024 * 8):
                        if chunk:
                            fp.write(chunk)
                            fp.flush()
                            os.fsync(fp.fileno())
            else:
                raise IOError(f"Failed to download SZ spectrum file: Error \
                                code {req.status_code}: {req.text}")

        sz_ell, sz_spec = np.loadtxt(self.sz_file, unpack=True, dtype=float)
        sz_ell = sz_ell.astype(int)

        # Making sure it's zero-indexed.
        self.sz_ell = np.arange(sz_ell.max() + 1)
        self.sz_spec = np.zeros_like(self.sz_ell)
        self.sz_spec[sz_ell] = sz_spec

        self.log.info("Loaded SZ.")

    def setup_for_tt_gibbs(self) -> None:
        if self.lowl_max < 2:
            raise LoggedError(self.log, f"lowl_max is set to {self.lowl_max}. \
                                          The Gibbs sampling requires \
                                          lowl_max >= 2.")

        gibbs_filename = os.path.normpath(
            os.path.join(self.data_folder, self.gibbs_sigma_filename))

        self.log.debug(f"Reading gibbs sigma data from from \
                         [{gibbs_filename}].")

        gibbs_file = fits.open(gibbs_filename)
        gibbs_chain = gibbs_file[0]

        lmax = gibbs_chain.header["LMAX"]
        numsamples = gibbs_chain.header["NUMSAMP"]
        numchains = gibbs_chain.header["NUMCHAIN"]
        numspec = gibbs_chain.header["NUMSPEC"]  # noqa F841

        data = gibbs_chain.data.T
        sigmas = data[:lmax, 0, :numchains, :numsamples]

        gibbs_file.close()

        self.log.debug("Done.")

        cl_filename = os.path.normpath(os.path.join(self.data_folder,
                                                    self.gibbs_cl_filename))

        self.log.debug(f"Reading gibbs cl data from from [{cl_filename}].")

        ls, cltt_fiducial = np.loadtxt(cl_filename, usecols=(0, 1),
                                       unpack=True)
        ls = ls.astype(int)
        self.cltt_fiducial = np.zeros((ls.max()+1, ))
        self.cltt_fiducial[ls] = cltt_fiducial

        self.log.debug("Done.")

        self.lmin = self.ttmin
        self.lmax = self.gibbs_ell_max

        self.sigmas = sigmas[:, :numchains,
                             self.gibbs_first_iteration-1:
                             self.gibbs_last_iteration:
                             self.gibbs_skip]
        self.numchain = self.sigmas.shape[1]
        self.numsamples = self.sigmas.shape[2]

        self.first_eval = True
        self.offset = -1.6375e30

        self.log.debug(f"Using values of ell from {self.lmin} to {self.lmax}.")
        self.log.debug(f"Using Gibbs samples from \
                        {self.gibbs_first_iteration} to \
                        {self.gibbs_last_iteration} in steps of \
                        {self.gibbs_skip}.")

        self.log.debug("Preparing Blackwell-Rao estimator using:")
        self.log.debug(f"num chains  = {self.numchain}")
        self.log.debug(f"num samples = {self.numsamples}")

        self.compute_br_estimator(self.cltt_fiducial)

        self.log.info("Initialized TT gibbs.")

    def setup_for_tt_beam_ptsrc(self, lmin: int, lmax: int) -> None:
        self.bptsrc_lmin = lmin
        self.bptsrc_lmax = lmax

        self.nmodes = 0
        self.ptsrc_mode_index = 0

        if self.beam_include_beam_modes:
            self.nmodes += self.n_beam_modes

        if self.beam_include_ptsrc_mode:
            self.ptsrc_mode_index = self.nmodes
            self.nmodes += 1

        if self.nmodes == 0:
            return

        self.beam_mode = np.zeros((lmax+1, self.n_beam_modes))
        self.fiducial_cltt = np.zeros((lmax+1, ))
        self.ptsrc_mode = np.zeros((lmax+1, ))

        if self.beam_include_beam_modes:
            ifn = os.path.normpath(os.path.join(self.data_folder,
                                                self.ifn_beam_modes))
            i, ls, x = np.loadtxt(ifn, usecols=(0, 1, 2), unpack=True)
            i = i.astype(int)
            ls = ls.astype(int)

            ls = ls[i <= self.n_beam_modes]
            x = x[i <= self.n_beam_modes]
            i = i[i <= self.n_beam_modes] - 1

            self.beam_mode[ls, i] = x

        if self.beam_fixed_fiducial_spectrum:
            ifn = os.path.normpath(os.path.join(self.data_folder,
                                                self.ifn_fiducial_cltt))
            ls, x = np.loadtxt(ifn, usecols=(0, 1), unpack=True)
            ls = ls.astype(int)

            self.fiducial_cltt[ls] = x

        if self.beam_include_ptsrc_mode:
            ifn = os.path.normpath(os.path.join(self.data_folder,
                                                self.ifn_ptsrc_mode))
            ls, x = np.loadtxt(ifn, usecols=(0, 1), unpack=True)
            ls = ls.astype(int)

            self.ptsrc_mode[ls] = x

        self.log.info("Initialized TT beam/ptsrc correction.")

    def setup_for_tt_exact(self) -> None:
        raise NotImplementedError("TT Pixel likelihood is not yet \
                                   implemented. Set .use_gibbs = True. or \
                                   disable the lowl TT likelihood.")

    def setup_for_tt_te_covariance(self) -> None:
        """
            Get TT/TE/TB diagonals
            Note that all the files are 2-indexed, but to make it easier on
            python coding (and to maintain the sanity of whomever needs to
            write/read this code), the datasets are re-ordered to be 0-indexed.
        """
        tt_filename = os.path.normpath(os.path.join(self.data_folder,
                                                    self.tt_filename))
        te_filename = os.path.normpath(os.path.join(self.data_folder,
                                                    self.te_filename))
        tb_filename = os.path.normpath(os.path.join(self.data_folder,
                                                    self.tb_filename))

        ttoff_filename = os.path.normpath(os.path.join(self.data_folder,
                                                       self.ttoff_filename))
        teoff_filename = os.path.normpath(os.path.join(self.data_folder,
                                                       self.teoff_filename))

        self.log.debug(f"Reading TT diagonal from [{tt_filename}].")

        ls, cltt_dat, ntt, fskytt = np.loadtxt(tt_filename, unpack=True,
                                               dtype=float)
        ls = ls.astype(int)
        self.cltt_dat = np.zeros((ls.max() + 1))
        self.ntt = np.zeros((ls.max() + 1))
        self.fskytt = np.zeros((ls.max() + 1))

        self.cltt_dat[ls] = cltt_dat
        self.ntt[ls] = ntt
        self.fskytt[ls] = fskytt

        self.log.debug(f"Read {self.cltt_dat.shape[0]} TT Cls.")
        self.log.debug(f"Read {self.ntt.shape[0]} TT n terms.")
        self.log.debug(f"Read {self.fskytt.shape[0]} TT sky terms.")

        self.log.debug(f"Reading TE diagonal from [{te_filename}].")

        ls, clte_dat, ntt_te, nee_te, fskyte = np.loadtxt(te_filename,
                                                          usecols=(0, 1, 3,
                                                                   4, 5),
                                                          unpack=True,
                                                          dtype=float)
        ls = ls.astype(int)
        self.clte_dat = np.zeros((ls.max() + 1))
        self.ntt_te = np.zeros((ls.max() + 1))
        self.nee_te = np.zeros((ls.max() + 1))
        self.fskyte = np.zeros((ls.max() + 1))

        self.clte_dat[ls] = clte_dat
        self.ntt_te[ls] = ntt_te
        self.nee_te[ls] = nee_te
        self.fskyte[ls] = fskyte

        self.log.debug(f"Read {self.clte_dat.shape[0]} TE Cls.")

        self.log.debug(f"Reading TB diagonal from [{tb_filename}].")

        ls, cltb_dat, ntt_tb, nbb_tb, fskytb = np.loadtxt(tb_filename,
                                                          usecols=(0, 1, 3,
                                                                   4, 5),
                                                          unpack=True,
                                                          dtype=float)
        ls = ls.astype(int)
        self.cltb_dat = np.zeros((ls.max() + 1))
        self.ntt_tb = np.zeros((ls.max() + 1))
        self.nbb_tb = np.zeros((ls.max() + 1))
        self.fskytb = np.zeros((ls.max() + 1))

        self.cltb_dat[ls] = cltb_dat
        self.ntt_tb[ls] = ntt_tb
        self.nbb_tb[ls] = nbb_tb
        self.fskytb[ls] = fskytb

        self.log.debug(f"Read {self.cltb_dat.shape[0]} TB Cls.")

        self.r_off_tttt = np.zeros((self.tt_lmax+1, self.tt_lmax+1))
        self.epsilon = np.zeros_like(self.r_off_tttt)
        self.r_off_tete = np.zeros((self.te_lmax+1, self.te_lmax+1))

        """ Get TT off-diagonal """
        i, j, epsilon, r_off_tttt = np.loadtxt(ttoff_filename,
                                               usecols=(0, 1, 2, 3),
                                               unpack=True)
        i = i.astype(int)
        j = j.astype(int)

        self.r_off_tttt[i, j] = r_off_tttt
        self.r_off_tttt[j, i] = r_off_tttt
        self.epsilon[i, j] = epsilon
        self.epsilon[j, i] = epsilon

        """ Get TE off-diagonal """
        i, j, r_off_tete = np.loadtxt(teoff_filename, usecols=(0, 1, 2),
                                      unpack=True)
        i = i.astype(int)
        j = j.astype(int)

        self.r_off_tete[i, j] = r_off_tete
        self.r_off_tete[j, i] = r_off_tete

        self.log.info("Initialized high-l TT/TE.")

    def load_fits_archive(self, filename):
        fp = fits.open(os.path.normpath(os.path.join(self.data_folder,
                                                     self.TEEEBB_maskfile)))
        a = np.asarray([x for x, _ in fp[1].data])
        b = np.asarray([y for _, y in fp[1].data])
        fp.close()
        return a, b

    def setup_for_lowl_te_ee_bb(self):
        """
        Initialize the TE/EE/BB low-l likelihood. Note that it is not
        implemented and should be disabled.
        """
        raise NotImplementedError("The Low-l polarization likelihood is not \
                                   implemented. Set .use_lowl_pol = False. to \
                                   disable this component. It may be \
                                   implemented in the future.")

        """
        Below here follows a start at copying the code from the
        original WMAP likelihood. We ended up not needing it (since there are
        much better low-l pol datasets available.
        """

        self.ires = 3
        self.nsmax = 2 ** self.ires
        self.nlmax = 3 * self.nsmax - 1
        self.np = 12 * self.nsmax ** 2

        self.mp = 0
        _, self.mask_r3 = self.load_fits_archive(self.TEEEBB_maskfile)

        self.ngood = np.zeros((self.np, ), dtype=int)
        for i in np.arange(self.np):
            if self.mask_r3[i] != 0:
                self.ngood[self.mp] = i
                self.mp += 1

        self.wl = np.loadtxt(
            os.path.normpath(os.path.join(self.data_folder,
                                          self.TEEEBB_filename_9)),
            usecols=(0, ))

        self.ninvplninv2 = np.zeros((self.mp * (2 * self.mp + 1),
                                     2 * self.nlmax))

        tmp = fits.open(os.path.normpath(os.path.join(self.data_folder,
                                                      "lowlP/std",
                                                      self.TEEEBB_filename_0)))
        ninvplninv3 = tmp[0].data.T
        tmp.close()

        k = 0
        for j in np.arange(self.mp * 2):
            for i in np.arange(j + 1):
                ip = self.ngood[i % self.mp] + self.np * (i // self.mp)
                jp = self.ngood[j % self.mp] + self.np * (j // self.mp)
                self.ninvplninv2[k, 1:self.nlmax-1] = \
                    ninvplninv3[2:self.nlmax, ip, jp]
                k += 1

        tmp = fits.open(os.path.normpath(os.path.join(self.data_folder,
                                                      "lowlP/std",
                                                      self.TEEEBB_filename_1)))
        ninvplninv3 = tmp[0].data.T
        tmp.close()

        k = 0
        for j in np.arange(self.mp * 2):
            for i in np.arange(j + 1):
                ip = self.ngood[i % self.mp] + self.np * (i // self.mp)
                jp = self.ngood[j % self.mp] + self.np * (j // self.mp)
                self.ninvplninv2[k, self.nlmax:2*self.nlmax-2] = \
                    ninvplninv3[2:self.nlmax, ip, jp]
                k += 1

        tmp = fits.open(os.path.normpath(os.path.join(self.data_folder,
                                                      "lowlP/std",
                                                      self.TEEEBB_filename_2)))
        NinvQUr3 = tmp[0].data.T
        tmp.close()

        self.w_r3 = np.zeros((2 * self.mp, ))

        ip = np.arange(self.mp)

        T, N = self.load_fits_archive(self.TEEEBB_filename_3)

        self.w_r3[ip] = T[self.ngood[ip]] * self.mask_r3[self.ngood[ip]]

        T, N = self.load_fits_archive(self.TEEEBB_filename_4)

        self.w_r3[ip + self.mp] = T[self.ngood[ip]] * \
            self.mask_r3[self.ngood[ip]]

        self.Dp0 = np.zeros((2 * self.mp, 2 * self.mp))

        ip = np.arange(self.mp)
        jp = np.arange(self.mp)
        ip, jp = np.meshgrid(ip, jp)

        self.Dp0[ip, jp] = NinvQUr3[self.ngood[ip], self.ngood[jp]]
        self.Dp0[ip, jp + self.mp] = \
            NinvQUr3[self.ngood[ip], self.ngood[jp] + self.np]
        self.Dp0[ip + self.mp, jp] = \
            NinvQUr3[self.ngood[ip] + self.np, self.ngood[jp]]
        self.Dp0[ip + self.mp, jp + self.mp] = \
            NinvQUr3[self.ngood[ip] + self.np, self.ngood[jp] + self.np]

        self.alm_tt = np.zeros((self.nlmax, self.nlmax))
        tmp_alm_tt = np.genfromtxt(
            os.path.normpath(os.path.join(self.data_folder,
                                          self.TEEEBB_filename_5)),
            converters={0: lambda x: np.complex(*np.safe_eval(x))})
        self.alm_tt[np.triu_indices(self.alm_tt.shape[0])] = \
            tmp_alm_tt[:self.nlmax * (self.nlmax + 1) // 2]

        tmp = fits.open(os.path.normpath(os.path.join(self.data_folder,
                                                      "lowlP/std",
                                                      self.TEEEBB_filename_6)))
        ninvy = tmp[0].data.T
        self.NinvY = np.zeros((ninvy.shape[1], ninvy.shape[2]), dtype=complex)
        self.NinvY[:, :] = ninvy[0, :, :] + 1.j * ninvy[1, :, :]
        tmp.close()

        self.log.info("Initialized low-l TE/EE/BB.")

    def setup_for_lowl_te_tb_ee_bb_eb(self):
        """ Initialize the TE/TB/EE/BB/EB low-l likelihood """
        raise NotImplementedError("TE/TB/EE/BB/EB low-l likelihood is not yet \
                                   implemented. Set .use_lowl_TBEB = False. \
                                   or .use_lowl_pol = False. to avoid this \
                                   error.")

    def compute_largest_term(self, cltt):
        self.log.info("First evaluation: Calculating largest term...")

        ls = np.arange(self.lmin, self.lmax+1)
        self.sigmas = self.sigmas[ls, ...]
        self.logsigma = np.log(self.sigmas)

        x = self.sigmas / cltt[ls, np.newaxis, np.newaxis]
        s = np.nansum(0.5 * (2.0 * ls[..., np.newaxis, np.newaxis] + 1) *
                      (np.log(x) - x) - self.logsigma, axis=0)

        self.offset = np.nanmax(s)

        self.log.info(f"Done finding largest term [{self.offset:g}].")

    def compute_br_estimator(self, cltt):
        """
          Compute the Blackwell-Rao estimator.
        """

        if self.first_eval:
            self.compute_largest_term(cltt)
            self.first_eval = False

        logl = 0.0
        ls = np.arange(self.lmin, self.lmax+1)

        x = self.sigmas / cltt[ls, np.newaxis, np.newaxis]
        s = np.nansum(0.5 * (2.0 * ls[..., np.newaxis, np.newaxis] + 1) *
                      (np.log(x) - x) - self.logsigma, axis=0)
        logl = np.log(np.nansum(np.exp(s - self.offset)))

        if logl < np.log(1e-20):
            logl = np.log(1e-30)

        return logl

    def get_requirements(self):
        return {"Cl": {k: 2000 for k in self.use_cl}}

    def loglike_lowl_TT(self, cltt):
        """
        Calculate the low l TT likelihood (with either the gibbs or
        pixel method).
        """

        logl = 0.0
        if self.use_gibbs:
            cltt_dummy = np.zeros((self.gibbs_ell_max+1, ))
            if self.gibbs_ell_max > self.lowl_max:
                cltt_dummy[self.lowl_max+1:self.gibbs_ell_max] = \
                    self.cltt_fiducial[self.lowl_max+1:self.gibbs_ell_max]

            cltt_dummy[self.ttmin:self.lowl_max+1] = \
                cltt[self.ttmin:self.lowl_max+1]

            logl = self.compute_br_estimator(cltt_dummy)
        else:
            raise NotImplementedError("Low l TT pixel likelihood not \
                                       implemented. Set .use_gibbs = True. or \
                                       disable low l TT likelihood.")

        return logl

    def loglike_lowl_TEEEBB(self, cltt, clte, clee, clbb):
        raise NotImplementedError("The Low-l polarization likelihood is not \
                                   implemented. Set .use_lowl_pol = False. to \
                                   disable this component. It may be \
                                   implemented in the future.")

        ls = 2 + np.arange(self.nlmax-1)

        lfac = (2.0 * np.pi) * 1e-6 * self.wl[ls] ** 2.0 / (ls * (ls + 1.0))

        xxx = np.zeros(2 * self.nlmax)

        xxx[ls-1] = (clee[ls] - clte[ls] ** 2.0 / cltt[ls]) * lfac
        xxx[ls-2 + self.nlmax] = clbb[ls] * lfac

        yyy = self.ninvplninv2 @ xxx
        Dp1 = np.zeros_like(self.Dp0)
        Dp1[np.tril_indices(Dp1.shape[0])] = yyy
        Dp1 = Dp1 + Dp1.T - np.diag(np.diag(Dp1))

        Dpu = linalg.cholesky(self.Dp0 + Dp1, lower=False)
        _, lndet = np.linalg.slogdet(Dpu)

        p_r3 = np.zeros((2 * self.mp, ))

        ip = np.arange(self.mp)
        ls = np.arange(2, 23)

        p_r3[ip] = np.nansum(clte[np.newaxis, ls] / cltt[np.newaxis, ls] *
                             self.wl[np.newaxis, ls] *
                             self.alm_tt[np.newaxis, ls, 0] *
                             self.NinvY[self.ngood[ip], 3, np.newaxis], axis=1)
        p_r3[ip + self.mp] = np.nansum(clte[np.newaxis, ls] /
                                       cltt[np.newaxis, ls] *
                                       self.wl[np.newaxis, ls] *
                                       self.alm_tt[np.newaxis, ls, 0] *
                                       self.NinvY[self.ngood[ip] + self.np,
                                                  3, np.newaxis],
                                       axis=1)

        chisqr = -0.5 * (self.w_r3 - p_r3) @ Dpu @ (self.w_r3 - p_r3)

        return chisqr, (self.TEEEBB_pixlike_logdet_offset - lndet * 2.0) / 2.0

    def loglike_lowl_TETBEEBBEB(self, cltt, clte, cltb, clee, clbb, cleb):
        return 0.0, 0.0

    def loglike_init_covterms(self, cls, TTTT, TTTE, TETE, TTTB, TBTB):
        """
        Calculate the TTTT/TTTE/TETE/TTTB/TBTB diagonal components of the
        high l TT/TE/TB likelihood.
        """
        ls = np.arange(self.ttmin, self.ttmax+1)

        if "tt" not in cls:
            cls["tt"] = np.zeros((self.ttmax+1))
        if "te" not in cls:
            cls["te"] = np.zeros((self.temax+1))
        if "tb" not in cls:
            cls["tb"] = np.zeros((self.temax+1))
        if "ee" not in cls:
            cls["ee"] = np.zeros((self.temax+1))
        if "bb" not in cls:
            cls["bb"] = np.zeros((self.temax+1))

        TTTT[ls] = 2.0 * np.power(cls["tt"][ls] + self.ntt[ls], 2.0) / \
            ((2.0 * ls + 1.0) * np.power(self.fskytt[ls], 2.0))

        ls = np.arange(self.temin, self.temax+1)

        TETE[ls] = ((cls["tt"][ls] + self.ntt_te[ls]) *
                    (cls["ee"][ls] + self.nee_te[ls]) +
                    np.power(cls["te"][ls], 2.0)) / \
            ((2.0 * ls + 1.0) * np.power(self.fskyte[ls], 2.0))
        TTTE[ls] = 2.0 * (cls["tt"][ls] + self.ntt[ls]) * cls["te"][ls] / \
            ((2.0 * ls + 1.0) * self.fskytt[ls] * self.fskyte[ls])

        TBTB[ls] = ((cls["tt"][ls] + self.ntt_tb[ls]) *
                    (cls["bb"][ls] + self.nbb_tb[ls]) +
                    np.power(cls["tb"][ls], 2.0)) / \
            ((2.0 * ls + 1.0) * np.power(self.fskytb[ls], 2.0))

    def get_fisher_TTTT(self, cltt, TTTT, TTTE, TETE):
        """
        Calculate the TTTT fisher matrix from the TTTT/TTTE/TETE components.
        """
        fisher = -self.r_off_tttt / np.sqrt(TTTT * TTTT[..., np.newaxis]) + \
            self.epsilon / (TTTT * TTTT[..., np.newaxis])

        ls = np.arange(self.ttmin, self.ttmax + 1)
        ls_te = np.arange(self.temin, self.te_lmax + 1)

        diag = np.zeros_like(TTTT)
        diag[ls] = 1.0 / TTTT[ls]
        diag[ls_te] = TETE[ls_te] / (TTTT[ls_te] * TETE[ls_te] -
                                     TTTE[ls_te] * TTTE[ls_te])
        np.fill_diagonal(fisher, diag)

        return fisher

    def loglike_highl_TT(self, cltt, TTTT, TTTE, TETE):
        off_log_curv = np.zeros((self.tt_lmax+1, self.tt_lmax+1))

        ls = np.arange(self.ttmin, self.ttmax + 1)
        z = np.log(self.cltt_dat[ls] + self.ntt[ls])
        zbar = np.log(cltt[ls] + self.ntt[ls])

        fisher = self.get_fisher_TTTT(cltt, TTTT, TTTE, TETE)

        # O[i, j] = ( cltt[i] + ntt[i] ) * fisher[i, j] * (cltt[j] + ntt[j])
        off_log_curv[self.ttmin:self.ttmax+1, self.ttmin:self.ttmax+1] = \
            (cltt[np.newaxis, ls] + self.ntt[np.newaxis, ls]) * \
            fisher[ls, :][:, ls] * \
            (cltt[ls, np.newaxis] + self.ntt[ls, np.newaxis])

        fisher[:self.tt_highl_start, :self.tt_highl_start] = 0.0
        off_log_curv[:self.tt_highl_start, :self.tt_highl_start] = 0.0

        dtt = cltt[ls] - self.cltt_dat[ls]

        logl = 0.0
        logl += (2.0 / 3.0) * \
            ((z - zbar) @ off_log_curv[ls, :][:, ls] @ (z - zbar))
        logl += (1.0 / 3.0) * (dtt @ fisher[ls, :][:, ls] @ dtt)

        return -logl / 2.0

    def loglike_TT_beam_ptsrc(self, cltt, fisher):
        """
        Add in the TT beam & point-source corrections term.
        """
        mode = np.zeros((self.bptsrc_lmax+1, self.nmodes))
        ls = np.arange(self.bptsrc_lmin, self.bptsrc_lmax + 1)

        if self.beam_include_beam_modes:
            if self.beam_fixed_fiducial_spectrum:
                mode[ls, :self.n_beam_modes] = self.beam_mode[ls, :] * \
                    self.fiducial_cltt[ls, np.newaxis]
            else:
                mode[ls, :self.n_beam_modes] = self.beam_mode[ls, :] * \
                    cltt[ls, np.newaxis]

        if self.beam_include_ptsrc_mode:
            mode[ls, self.ptsrc_mode_index] = self.ptsrc_err * \
                self.ptsrc_mode[ls]

        a = np.zeros((self.nmodes, ))
        b = np.identity(self.nmodes)
        F_mode = np.zeros((self.bptsrc_lmax+1, self.nmodes))

        if self.beam_diagonal_sigma:
            F_mode[ls, :] = np.diag(fisher)[ls, np.newaxis] * mode[ls, :]
        else:
            # I havent found a good way to unloop this.
            for l1 in ls:
                for l2 in ls:
                    if l2 < l1:
                        F_mode[l1, :] += fisher[l2, l1] * mode[l2, :]
                    else:
                        F_mode[l1, :] += fisher[l1, l2] * mode[l2, :]

        a[:self.nmodes] = np.nansum(
            (self.cltt_dat[ls] - cltt[ls])[..., np.newaxis] * F_mode[ls, :],
            axis=0)
        b += np.dot(mode.T, F_mode)

        b = linalg.cholesky(b, lower=True)
        c = linalg.solve_triangular(b, a, lower=True)
        dgauss = np.dot(a, c)

        if self.beam_gaussian_likelihood:
            dlndet = np.log(np.nanprod(np.diag(b) ** 2.0))

            return (dgauss - dlndet) / 2.0

        return -dgauss

    def get_fisher_TETE(self, clte, TTTT, TTTE, TETE):
        """
        Calculate the TETE fisher matrix from the TTTT/TTTE/TETE components.
        """
        fisher = -self.r_off_tete / np.sqrt(TETE * TETE[..., np.newaxis])

        ls = np.arange(self.temin, self.temax + 1)
        ls_lo = ls[ls <= self.te_lmax]

        diag = np.zeros((fisher.shape[0]))
        diag[ls] = 1.0 / TETE[ls]
        diag[ls_lo] = TTTT[ls_lo] / (TTTT[ls_lo] * TETE[ls_lo] -
                                     TTTE[ls_lo] * TTTE[ls_lo])
        np.fill_diagonal(fisher, diag)

        return fisher

    def loglike_highl_TE(self, clte, TTTT, TTTE, TETE):
        """
        Calculate the TETE log-likelihood.
        Returns both the model chi2 and the log(det) of the fisher matrix.
        """
        ls = np.arange(self.temin, self.temax + 1)

        fisher = self.get_fisher_TETE(clte, TTTT, TTTE, TETE)
        fisher[:self.te_highl_start, :self.te_highl_start] = 0.0

        dte = clte[ls] - self.clte_dat[ls]

        logl = (dte @ fisher[ls, :][:, ls] @ dte)
        _, logdet = np.linalg.slogdet(fisher[self.te_highl_start:,
                                             self.te_highl_start:])

        return -logl / 2.0, (logdet + self.TE_logdet_offset) / 2.0

    def loglike_highl_TB(self, cltb, TTTT, TTTB, TBTB):
        """
        Calculate the TBTB log-likelihood.
        Returns both the model chi2 and the log(det) of the fisher matrix.
        """
        ls = np.arange(self.temin, self.temax + 1)

        fisher = np.zeros((self.te_lmax+1, self.te_lmax+1))

        diag = np.zeros((fisher.shape[0]))
        diag[ls] = 1.0 / TBTB[ls]
        np.fill_diagonal(fisher, diag)
        fisher[:self.te_highl_start, :self.te_highl_start] = 0.0

        dtb = cltb[ls] - self.cltb_dat[ls]

        logl = (dtb @ fisher[ls, :][:, ls] @ dtb)
        _, logdet = np.linalg.slogdet(fisher[self.te_highl_start:,
                                             self.te_highl_start:])

        return -logl / 2.0, (logdet + self.TB_logdet_offset) / 2.0

    def loglike(self, Cls, **params):
        components = {}

        if self.use_sz:
            sz_spec = np.zeros((Cls.get("ell").max()+1, ))
            lmax = min(len(sz_spec), len(self.sz_spec))
            sz_spec[:lmax] = self.sz_spec[:lmax]
            Cls["tt"] = Cls.get("tt") + params["A_sz"] * \
                sz_spec[Cls.get("ell")]

        if self.use_lowl_TT:
            components["lowl_TT_gibbs"] = self.loglike_lowl_TT(Cls.get("tt"))

        TTTT = np.zeros((self.tt_lmax+1, ))
        TETE = np.zeros((self.te_lmax+1, ))
        TTTE = np.zeros((self.te_lmax+1, ))
        TBTB = np.zeros((self.te_lmax+1, ))
        TTTB = np.zeros((self.te_lmax+1, ))

        if self.use_highl_TT or self.use_highl_TE or self.use_highl_TB:
            self.loglike_init_covterms(Cls, TTTT, TTTE, TETE, TTTB, TBTB)

        if self.use_highl_TT:
            components["MASTER_TTTT"] = self.loglike_highl_TT(Cls.get("tt"),
                                                              TTTT, TTTE, TETE)

        if self.use_highl_TE:
            logchi, logdet = self.loglike_highl_TE(Cls.get("te"), TTTT, TTTE,
                                                   TETE)
            components["MASTER_TETE_chi2"] = logchi
            components["MASTER_TETE_det"] = logdet

        if self.use_highl_TB:
            logchi, logdet = self.loglike_highl_TB(Cls.get("tb"), TTTT, TTTB,
                                                   TBTB)
            components["MASTER_TBTB_chi2"] = logchi
            components["MASTER_TBTB_det"] = logdet

        if self.use_highl_TT and self.use_highl_TT_beam_ptsrc:
            fisher = self.get_fisher_TTTT(Cls.get("tt"), TTTT, TTTE, TETE)
            components["beamptsrc_TT"] = self.loglike_TT_beam_ptsrc(
                Cls.get("tt"), fisher)

        logp = 0.0
        for c in components:
            logp += components[c]

        return logp, components

    def calculate(self, state, want_derived=True, **params):
        Cls = self.provider.get_Cl(ell_factor=True)

        logp, components = self.loglike(Cls, **params)

        state["logp"] = logp
        if state["derived"] is None:
            state["derived"] = {}
        for c in components:
            state["derived"]["chi2_wmap_" + c] = -2.0 * components[c]
            state["derived"]["logp_" + c] = components[c]

        chi2 = -2.0 * logp

        self.log.debug(f"logL = {logp:g}, Chisqr = {chi2:.2f}")
        for c in components:
            lp = components[c]
            chi2 = -2.0 * lp
            self.log.debug(f"\t{c}: {lp:g}, Chisqr = {chi2:.2f}")

        return True

    def logp(self, **params):
        Cls = self.provider.get_Cl(ell_factor=True)
        loglike, _ = self.loglike(Cls, **params)

        return loglike


class WMAPBestFitv5(Theory):
    def initialize(self):
        data_file_path = os.path.normpath(getattr(self, "path", None) or
                                          os.path.join(self.packages_path,
                                                       "data"))
        self.data_folder = os.path.join(data_file_path, self.data_folder)

        if not os.path.exists(self.data_folder):
            raise LoggedError(self.log, f"No data folder found at \
                                          [{self.data_folder}].")

        self.cl_file = os.path.join(self.data_folder, self.cl_filename)

        cldata = np.loadtxt(self.cl_file)

        self.cls_dict = {
            k: np.zeros((cldata.shape[0]+2, )) for k in ["ell", "tt", "te",
                                                         "ee", "tb", "eb",
                                                         "bb"]}

        ls = cldata[:, 0].astype(int)
        self.cls_dict["ell"] = np.arange(
            self.cls_dict["ell"].shape[0]).astype(int)
        self.cls_dict["tt"][ls] = cldata[:, 1]
        self.cls_dict["ee"][ls] = cldata[:, 2]
        self.cls_dict["bb"][ls] = cldata[:, 3]
        self.cls_dict["te"][ls] = cldata[:, 4]

        self.log.info("Initialized WMAP best fit v5 model.")

    def get_can_provide(self):
        return ["Cl"]

    def get_Cl(self, ell_factor=True, units="1"):
        ls_fac = np.ones_like(self.cls_dict["ell"])

        if not ell_factor:
            ls_fac = (2.0 * np.pi) / \
                (self.cls_dict["ell"] * (self.cls_dict["ell"] + 1.0))

        cmb_fac = self._cmb_unit_factor(units, 2.726)

        cls = {
            k: self.current_state["cl_amp"] * self.cls_dict[k] * ls_fac *
            cmb_fac ** 2.0 for k in self.cls_dict
        }

        cls["ell"] = self.cls_dict["ell"]

        return cls

    def _cmb_unit_factor(self, units: str, T_cmb: float) -> float:
        units_factors = {
            "1": 1,
            "muK2": T_cmb * 1.0e6,
            "K2": T_cmb,
            "FIRASmuK2": 2.7255e6,
            "FIRASK2": 2.7255
        }

        try:
            return units_factors[units]
        except KeyError:
            raise LoggedError(self.log, f"Units '{units}' not recognized. Use \
                                          one of {list(units_factors)}.")

    def calculate(self, state: dict, want_derived: bool = True,
                  **params) -> None:
        state["cl_amp"] = params["cl_amp"]

    def get_can_support_params(self) -> Sequence[str]:
        return ["cl_amp"]
