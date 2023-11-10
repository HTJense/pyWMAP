import pytest  # noqa F401
import numpy as np
from cobaya.model import get_model

info = {
    "params": {
        "cl_amp": 1.0
    },
    "theory": {
        "wmaplike.WMAPBestFitv5": None
    },
    "likelihood": {
        "wmaplike.WMAPLike": {
            "use_lowl_TT": True,
            "use_highl_TT": True,
            "use_highl_TE": True,
            "use_highl_TB": False,
            "use_lowl_pol": False,
            "use_lowl_TBEB": False,
            "use_highl_TT_beam_ptsrc": False,
            "use_sz": False,
        }
    },
    "sampler": {"evaluate": None},
    "debug": False
}

chisqr_expected = {
    "lowl_TT_gibbs": -13.614869,
    "MASTER_TTTT": 1200.005224,
    "MASTER_TETE_chi2": 831.787122,
    "MASTER_TETE_det": 3605.526233,
    "MASTER_TBTB_chi2": 775.110063,
    "MASTER_TBTB_det": 3610.385517,
}


def test_import():
    import wmaplike  # noqa F401


def test_model():
    model = get_model(info)  # noqa F841


def test_chi2():
    info["likelihood"] = {
        "wmaplike.WMAPLike": {
            "use_lowl_TT": True,
            "use_highl_TT": True,
            "use_highl_TE": True,
            "use_highl_TB": True,
            "use_lowl_pol": False,
            "use_lowl_TBEB": False,
            "use_highl_TT_beam_ptsrc": False,
            "use_sz": False,
        }
    }

    model = get_model(info)
    model.logposterior({"cl_amp": 1.0, "A_sz": 0.0})
    like = model.likelihood["wmaplike.WMAPLike"]
    chi2s = like.current_state["derived"]

    for component in chisqr_expected:
        if f"chi2_wmap_{component}" in chi2s:
            assert np.isclose(chi2s[f"chi2_wmap_{component}"],
                              chisqr_expected[component])
            print(f"{component} is fine.")


if __name__ == "__main__":
    test_import()
    test_model()
    test_chi2()
