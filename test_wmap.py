import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from wmaplike import WMAPLike

packages_path = os.environ.get("COBAYA_PACKAGES_PATH")
if packages_path is None:
	packages_path = "/home/ojima/Cardiff/cobaya/packages"

like = WMAPLike({
	"packages_path" : packages_path,
	"use_lowl_TT" : True,
	"use_highl_TT" : True,
	"use_highl_TE" : True,
	"use_highl_TB" : False,
	"use_lowl_pol" : False,
	"use_lowl_TBEB" : False,
	"debug_timing" : True
})

bestfit = np.loadtxt(os.path.normpath(os.path.join(packages_path, "data/WMAP/dr5/wmap_likelihood_v5/data/test_cls_v5.dat")))

cls_dict = { k : np.zeros((bestfit.shape[0]+2,)) for k in [ "tt", "te", "ee", "tb", "eb", "bb" ] }

ls = bestfit[:,0].astype(int)
cls_dict["tt"][ls] = bestfit[:,1]
cls_dict["ee"][ls] = bestfit[:,2]
cls_dict["bb"][ls] = bestfit[:,3]
cls_dict["te"][ls] = bestfit[:,4]

print("‒" * 80)
print("\tTesting gibbs TT code")

offset_exp = -759.901445147943

print(f"\tCalculated offset is {like.offset:.5f}")
print(f"\tExpected offset is {offset_exp:.5f}")

#sys.exit()

print("‒" * 80)
print("\tTesting MASTER TTTT code")

TTTT = np.zeros((like.tt_lmax+1,))
TETE = np.zeros((like.te_lmax+1,))
TTTE = np.zeros((like.te_lmax+1,))
TTTB = np.zeros((like.te_lmax+1,))
TBTB = np.zeros((like.te_lmax+1,))

like.loglike_init_covterms(cls_dict, TTTT, TTTE, TETE, TTTB, TBTB)

if False:
	fig, ax = plt.subplots(2, 2, figsize = (12, 8))

	ls, TTTT_baseline = np.loadtxt("data/tttt", dtype = float, unpack = True)
	ls = ls.astype(int)

	ax[0,0].plot(ls, TTTT[ls], lw = 2, c = "C0")
	ax[0,0].plot(ls, TTTT_baseline, lw = 2, ls = ":", c = "k")
	ax[0,0].semilogy()

	ls, TTTE_baseline, TETE_baseline, TBTB_baseline = np.loadtxt("data/ttte_tete_tbtb", dtype = float, unpack = True)
	ls = ls.astype(int)

	ax[0,1].plot(ls, TTTE[ls], lw = 2, c = "C1")
	ax[0,1].plot(ls, TTTE_baseline, lw = 2, ls = ":", c = "k")

	ax[1,0].plot(ls, TETE[ls], lw = 2, c = "C2")
	ax[1,0].plot(ls, TETE_baseline, lw = 2, ls = ":", c = "k")
	ax[1,0].semilogy()

	ax[1,1].plot(ls, TBTB[ls], lw = 2, c = "C3")
	ax[1,1].plot(ls, TBTB_baseline, lw = 2, ls = ":", c = "k")

	plt.show()

fisher = like.get_fisher_TTTT(cls_dict.get("tt"), TTTT, TTTE, TETE)

fisher_baseline = np.loadtxt("data/fisher_tttt", dtype = float)

# triangular into symmetrical matrix.
fisher_baseline = fisher_baseline + fisher_baseline.T - np.diagflat(np.diag(fisher_baseline))

print("‒" * 80)

avg_diff = np.nanmean( np.abs(fisher[2:,2:] - fisher_baseline) )
max_diff = np.nanmax( np.abs(fisher[2:,2:] - fisher_baseline) )
max_diff_idx = np.abs(fisher[2:,2:] - fisher_baseline).argmax()
max_diff_val = fisher_baseline.flatten()[max_diff_idx]

print(f"\tTTTT MASTER FISHER")
print(f"\tMean absolute difference: {avg_diff:.3e}.")
print(f"\tMaximum absolute difference: {max_diff:.3e} [{max_diff_val:.4e}].")
print(f"\tDifferences of 10^-10 or less are to be expected due to floating-point precision in the baseline file.")

if False:
	fig, ax = plt.subplots(2, 3, figsize = (12, 8), sharex = True, sharey = True)

	im = ax[0,0].imshow( fisher[2:,2:] )
	fig.colorbar(im, ax = ax[0,0])
	im = ax[0,1].imshow( fisher_baseline )
	fig.colorbar(im, ax = ax[0,1])
	im = ax[0,2].imshow( fisher[2:,2:] - fisher_baseline )
	fig.colorbar(im, ax = ax[0,2])

	im = ax[1,0].imshow( np.log10(np.abs(fisher[2:,2:])) )
	fig.colorbar(im, ax = ax[1,0])
	im = ax[1,1].imshow( np.log10(np.abs(fisher_baseline)) )
	fig.colorbar(im, ax = ax[1,1])
	im = ax[1,2].imshow( np.log10(np.abs( fisher[2:,2:] - fisher_baseline )) )
	fig.colorbar(im, ax = ax[1,2])

	plt.show()

ls, z_baseline, zbar_baseline, zdiff_baseline = np.loadtxt("data/ttvecs", usecols = (0,1,2,3), unpack = True, dtype = float)

ls = ls.astype(int)
z = np.log(like.cltt_dat[ls] + like.ntt[ls])
zbar = np.log(cls_dict["tt"][ls] + like.ntt[ls])
zdiff = z - zbar

z_diff = np.nanmean( np.log10(np.abs(z_baseline - z)) )
zbar_diff = np.nanmean( np.log10(np.abs(zbar_baseline - zbar)) )
zdiff_diff = np.nanmean( np.log10(np.abs(zdiff_baseline - zdiff)) )

print("")
print(f"\tTTTT MASTER VECTORS")
print(f"\tSum log difference z    : {z_diff:.3e}.")
print(f"\tSum log difference zbar : {zbar_diff:.3e}.")
print(f"\tSum log difference zdiff: {zdiff_diff:.3e}.")

if False:
	fig, ax = plt.subplots(2, 3, figsize = (12, 8))
	
	ax[0,0].plot(ls, z)
	ax[0,0].plot(ls, z_baseline, c = "k", ls = ":")
	
	ax[0,1].plot(ls, zbar)
	ax[0,1].plot(ls, zbar_baseline, c = "k", ls = ":")
	
	ax[0,2].plot(ls, zdiff)
	ax[0,2].plot(ls, zdiff_baseline, c = "k", ls = ":")
	
	ax[1,0].plot(ls, z - z_baseline)
	ax[1,0].plot(ls, z_baseline - z_baseline, c = "k", ls = ":")
	
	ax[1,1].plot(ls, zbar - zbar_baseline)
	ax[1,1].plot(ls, zbar_baseline - zbar_baseline, c = "k", ls = ":")
	
	ax[1,2].plot(ls, zdiff - zdiff_baseline)
	ax[1,2].plot(ls, zdiff_baseline - zdiff_baseline, c = "k", ls = ":")
	
	plt.show()

ls, cltt_temp_baseline, cltt_dat_baseline, cltt_diff_baseline = np.loadtxt("data/ttvecs", usecols = (0,4,5,6), unpack = True, dtype = float)
ls = ls.astype(int)

cltt_temp = cls_dict["tt"][ls]
cltt_dat = like.cltt_dat[ls]
cltt_diff = cltt_temp - cltt_dat

if False:
	fig, ax = plt.subplots(2, 3, figsize = (12, 8))
	
	ax[0,0].plot(ls, cltt_temp)
	ax[0,0].plot(ls, cltt_temp_baseline, c = "k", ls = ":")
	
	ax[0,1].plot(ls, cltt_dat)
	ax[0,1].plot(ls, cltt_dat_baseline, c = "k", ls = ":")
	
	ax[0,2].plot(ls, cltt_diff)
	ax[0,2].plot(ls, cltt_diff_baseline, c = "k", ls = ":")
	
	ax[1,0].plot(ls, cltt_temp - cltt_temp_baseline)
	ax[1,0].plot(ls, cltt_temp_baseline - cltt_temp_baseline, c = "k", ls = ":")
	
	ax[1,1].plot(ls, cltt_dat - cltt_dat_baseline)
	ax[1,1].plot(ls, cltt_dat_baseline - cltt_dat_baseline, c = "k", ls = ":")
	
	ax[1,2].plot(ls, cltt_diff - cltt_diff_baseline)
	ax[1,2].plot(ls, cltt_diff_baseline - cltt_diff_baseline, c = "k", ls = ":")
	
	plt.show()

print("‒" * 80)

print(f"\tTETE MASTER FISHER")

fisher = like.get_fisher_TETE(cls_dict.get("te"), TTTT, TTTE, TETE)

fisher_baseline = np.loadtxt("data/fisher_tete", dtype = float)
# triangular into symmetrical matrix.
fisher_baseline = fisher_baseline + fisher_baseline.T - np.diagflat(np.diag(fisher_baseline))

avg_diff = np.nanmean( np.abs(fisher[2:,2:] - fisher_baseline) )
max_diff = np.nanmax( np.abs(fisher[2:,2:] - fisher_baseline) )
max_diff_idx = np.abs(fisher[2:,2:] - fisher_baseline).argmax()
max_diff_val = fisher_baseline.flatten()[max_diff_idx]

print(f"\tMean absolute difference: {avg_diff:.3e}.")
print(f"\tMaximum absolute difference: {max_diff:.3e} [{max_diff_val:.4e}].")
print(f"\tDifferences of 10^-8  or less are to be expected due to floating-point precision in the baseline file.")

if False:
	fig, ax = plt.subplots(2, 3, figsize = (12, 8), sharex = True, sharey = True)
	
	im = ax[0,0].imshow( fisher[2:,2:] )
	fig.colorbar(im, ax = ax[0,0])
	
	im = ax[0,1].imshow( np.log(np.abs(fisher[2:,2:])) )
	fig.colorbar(im, ax = ax[0,1])
	
	im = ax[0,2].imshow( fisher[2:,2:] - fisher_baseline )
	fig.colorbar(im, ax = ax[0,2])
	
	im = ax[1,0].imshow( fisher_baseline )
	fig.colorbar(im, ax = ax[1,0])
	
	im = ax[1,1].imshow( np.log(np.abs(fisher_baseline)) )
	fig.colorbar(im, ax = ax[1,1])
	
	im = ax[1,2].imshow( np.log10(np.abs( fisher[2:,2:] - fisher_baseline )) )
	fig.colorbar(im, ax = ax[1,2])
	
	plt.show()

logp_tot, components = like.loglike(cls_dict)

chisqr_expected = {
	"lowl_TT_gibbs" : -13.614869,
	"MASTER_TTTT" : 1200.005224,
	"MASTER_TETE_chi2" : 815.433752,
	"MASTER_TETE_det" : 3541.537184,
	"MASTER_TBTB_chi2" : 756.521494,
	"MASTER_TBTB_det" : 3542.432360,
	"TEEEBB_chi2" : 1320.994614,
	"TEEEBB_det" : 692.874562,
	"beamptsrc_TT" : 0.735565
}

#chisqr_tot_expected = 7557.965820
logp_tot_avail = sum([ components[c] if c in chisqr_expected else 0.0 for c in components ])
chisqr_tot_expected = sum([ chisqr_expected[c] if c in chisqr_expected else 0.0 for c in components ])

is_fine = { k : False for k in components }

print("‒" * 80)
print("Breakdown of loglikes")
print("‒" * 80)
print("Parameters")
print("ttmax - lowl_max = " + str(like.ttmax - like.lowl_max))
if like.use_gibbs: print(f"lmin : lmax = {like.lmin} : {like.lmax}")
print("‒" * 80)
for c in components:
	logp = components[c]
	chi2 = -2.0 * logp
	chi2_exp = chisqr_expected[c] if c in chisqr_expected else None
	print(f"{c:19s}: logL = {logp:8g}, chisqr = {chi2:12.6f}" + (f", exp = {chi2_exp:12.6f}" if chi2_exp is not None else ""))
	if chi2_exp is not None and np.abs(chi2 - chi2_exp) >= 1.0:
		print(f"\tMajor differences detected!")
	elif chi2_exp is not None and np.abs(chi2 - chi2_exp) >= 0.01:
		print(f"\tMinor differences detected.")
	else:
		is_fine[c] = True
print("‒" * 80)

chi2_tot = -2.0 * logp_tot
chi2_tot_avail = -2.0 * logp_tot_avail

print(f"TOTAL: logL = {logp_tot:g}, chisqr = {chi2_tot:.3f}")
print(f"COMPARISON: logL = {logp_tot_avail:g}, chisqr = {chi2_tot_avail:.3f}, exp = {chisqr_tot_expected:.3f}")

print("‒" * 80)

found_errors = False

for c in components:
	if not is_fine[c]:
		found_errors = True
		print(f"\tCheck component {c}.")

if not found_errors:
	print(f"\tNo major problems detected! Enjoy WMAP.")

print("‒" * 80)
