#!/usr/bin/env python
# coding: utf-8

# In[2]:


# ── config.py  (Python ≤3.9 compatible) ──────────────────────────────
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Sequence, Any, Optional, Dict, List


# ───────────────── Spike detection parameters ────────────────────────
@dataclass
class SpikeConfig:
    detmethod: Any                         # e.g. AP_detection_lofmad
    fs: int = 20_000
    cutoff: int = 100
    order: int = 9
    dead_time: float = 0.5
    threshold: float = 5.5
    usewelllist: Optional[Sequence[int]] = None
    negative: bool = True
    positive: bool = True
    PipNoise: bool = False
    noallwells: bool = False
    mindur: int = 3
    maxdur: int = 500
    E_Stim: bool = False
    MetaMaskPath: Optional[str] = None

    # Convenience helper if a legacy function still needs **kwargs
    def to_kwargs(self) -> dict:
        d = asdict(self)
        d.pop("detmethod")
        return d


# ───────────────── Burst detection parameters ────────────────────────
@dataclass
class BurstConfig:
    """
    Parameters for burst detection.

    Attributes
    ----------
    method : str
        One of 'maxinterval', 'logisi', 'poisson'.
    param_set : str
        Human-readable label, e.g. 'set7'.
    values : List[int]
        Raw numbers that your detector expects.
        * maxinterval → [startISI, endISI, minIBI, minNSpikes, …]
        * logisi      → [threshold_us]
        * poisson     → [...]
    """
    method: str = "maxinterval"
    param_set: str = "default"
    values: List[int] = None          # list length depends on method

    def to_kwargs(self) -> Dict[str, Dict[str, List[int]]]:
        """
        Convert to the legacy kwargs shape:
        {'maxinterval': {'set7': [...]}}
        """
        return {self.method: {self.param_set: self.values}}


# In[ ]:




