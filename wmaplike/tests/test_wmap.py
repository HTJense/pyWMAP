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
            "debug": True
        }
    },
    "sampler": {"evaluate": None},
    "debug": True
}

chisqr_expected = {
    "lowl_TT_gibbs" : -13.614869,
    "MASTER_TTTT" : 1200.005224,
    "MASTER_TETE_chi2": 831.787122,
    "MASTER_TETE_det": 3605.526233,
    "MASTER_TBTB_chi2": 775.110063,
    "MASTER_TBTB_det": 3610.385517,
    "TEEEBB_chi2" : 1320.994614,
    "TEEEBB_det" : 692.874562,
    "beamptsrc_TT" : 0.735565
}

def test_import():
    import wmaplike  # noqa F401


def test_model():
    model = get_model(info)  # noqa F841


def test_highl():
    info["likelihood"] = {
        "wmaplike.WMAPLike": {
            "use_lowl_TT": True,
            "use_highl_TT": True,
            "use_highl_TE": True,
            "use_highl_TB": False,
            "use_lowl_pol": False,
            "use_lowl_TBEB": False,
            "use_highl_TT_beam_ptsrc": False,
            "use_sz": False,
            "debug": True
        }
    }

    model = get_model(info)
    model.loglikes()
    
