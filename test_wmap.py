import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from wmaplike import WMAPLike

packages_path = os.environ.get("COBAYA_PACKAGES_PATH")
if packages_path is None:
	raise Exception("Cannot locate Cobaya packages path. Are you sure it is set?")

like = WMAPLike({
	"packages_path" : packages_path,
	"use_lowl_TT" : True,
	"use_highl_TT" : True,
	"use_highl_TE" : True,
	"use_highl_TB" : False,
	"use_lowl_pol" : False,
	"use_lowl_TBEB" : False
})

bestfit = np.loadtxt(os.path.normpath(os.path.join(packages_path, like.data_folder, "data/test_cls_v5.dat")))

cls_dict = { k : np.zeros((bestfit.shape[0]+2,)) for k in [ "tt", "te", "ee", "tb", "eb", "bb" ] }

ls = bestfit[:,0].astype(int)
cls_dict["tt"][ls] = bestfit[:,1]
cls_dict["ee"][ls] = bestfit[:,2]
cls_dict["bb"][ls] = bestfit[:,3]
cls_dict["te"][ls] = bestfit[:,4]

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
