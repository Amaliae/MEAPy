#!/usr/bin/env python
# coding: utf-8

# In[6]:


"""
src/meafea/spike_detection.py
───────────────────────────────────────────────────────────────────────
Self-contained MEA spike-detection toolbox + serial folder runner.
Python 3.8 compatible.
"""

# --------------------------------------------------------------------
# Standard / third-party imports
# --------------------------------------------------------------------
from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional, Sequence
from dataclasses import dataclass, asdict

import os, re

import numpy as np
import pandas as pd
from lxml import etree
import McsPy.McsData
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation, skew, kurtosis
from sklearn.decomposition import PCA
from datetime import datetime as dt 
import meafea.moduleV15_vs1 as m15 


# --------------------------------------------------------------------
# 5.  Parameter bundle for serial runner
# --------------------------------------------------------------------
@dataclass
class SpikeConfig:
    # --- 1.  positional parameters (legacy order) -------------------
    identifier1: str = "h5"                 # substring inside filename
    identifier2: str = "h5"                 #   (kept for API compat)
    identifier3: str = "downsampled_data"   # folder / file token
    MetadataMaskPath: str = ""
    threshold_dict: Dict[str, Any] | None = None   # {'File': ..., 'Dose Label': ...}
    Analyze: str = "True"                   # "True" | "False"

    # --- 2.  keyword parameters -------------------------------------
    detmethod: Any = None                  # your detection function
    onlyLed: bool = True
    E_Stim: bool = False
    fs: int = 20_000
    cutoff: int = 100
    highpassth: int = 5
    order: int = 9
    dead_time: float = 0.5
    threshold: float = 4.5
    usewelllist: Optional[Sequence[int]] = None
    negative: bool = True
    positive: bool = True
    PipNoise: bool = False
    noallwells: bool = False
    mindur: int = 3
    maxdur: int = 500
    detsortfeatures: str = "False"

    # ----------------------------------------------------------------
    # helpers
    # ----------------------------------------------------------------
    def pos_args(self) -> List[Any]:
        """
        Return the six positional arguments
            (identifier1, identifier2, identifier3,
             MetadataMaskPath, threshold_dict, Analyze)
        in the exact order that your legacy
        experiment_spike_detection(...) expects.
        """
        return [
            self.identifier1,
            self.identifier2,
            self.identifier3,
            self.MetadataMaskPath,
            self.threshold_dict or {},
            self.Analyze,
        ]

    def kw_args(self) -> Dict[str, Any]:
        """
        All remaining parameters as **kwargs** (baseline for
        experiment_spike_detection).
        """
        d = asdict(self)
        for k in (
            "identifier1", "identifier2", "identifier3",
            "MetadataMaskPath", "threshold_dict", "Analyze"
        ):
            d.pop(k)
        return d

    # (Optional) convenience – full dict view
    def asdict(self) -> Dict[str, Any]:
        return asdict(self)
# -------------------------------------------------------------------
# helper: SpikeConfig  ->  (positional_args, keyword_args)
# -------------------------------------------------------------------
def _config_to_positional(cfg: "SpikeConfig"):
    """
    Return
        pos_args : list
            [identifier1, identifier2, identifier3,
             MetadataMaskPath, threshold_dict, Analyze]
        kw_args  : dict
            every remaining field of SpikeConfig
    """
    pos_args = cfg.pos_args()      # uses the method we just added
    kw_args  = cfg.kw_args()       # the remainder
    return pos_args, kw_args


def run_spike_detection_serial(
    hdf5_root : str | Path,
    save_root : str | Path,
    cfg       : "SpikeConfig",
) -> None:
    """
    Call *experiment_spike_detection* ONCE, passing the entire folder.

    Parameters
    ----------
    hdf5_root : directory that contains all *.h5 / *.mwd recordings
    save_root : directory where <file>Spikes.csv and stamp.npy will go
    cfg       : SpikeConfig (tunable parameters)
    """
    from meafea.moduleV15_vs1 import experiment_spike_detection

    hdf5_root = Path(hdf5_root).expanduser().resolve()
    save_root = Path(save_root).expanduser().resolve()
    save_root.mkdir(parents=True, exist_ok=True)

    # Convert config → positional / keyword bundles
    pos_args, kw_args = _config_to_positional(cfg)

    print(f"[{dt.now():%H:%M:%S}]  spike detection on folder  {hdf5_root}")

    experiment_spike_detection(
       os.path.dirname(hdf5_root),   # <-- whole directory
       save_root,   # where to write *Spikes.csv
        *pos_args,               # identifier1, identifier2, ...
        **kw_args,               # detmethod, fs, threshold, ...
    )

    print("✓ spike detection finished →", save_root)
