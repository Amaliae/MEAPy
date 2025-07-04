#!/usr/bin/env python
# coding: utf-8

# In[2]:


from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Sequence, Union
import pandas as pd
import numpy as np
from datetime import datetime

@dataclass
class BurstConfig:
    detection_methods : List[str]
    parameter_dict    : Dict[str, Dict[str, List[int]]]

    # convenience helpers
    def args(self):   return [self.detection_methods]
    def kwargs(self): return self.parameter_dict


def run_burst_detection_serial(
    spike_root : Union[str, Path],
    save_root  : Union[str, Path],
    cfg        : "BurstConfig",
) -> None:
    """
    Calls your original experiment_burst_detection once, passing
        (folder_path, save_to, identifier1, identifier2, method_set, params_set)
    exactly as defined in BurstDetectionPossion_lib_vs1.
    """
    from meafea.BurstDetection_lib_vs1 import experiment_burst_detection

    spike_root = Path(spike_root).expanduser().resolve()
    save_root  = Path(save_root).expanduser().resolve()
    save_root.mkdir(parents=True, exist_ok=True)
    

    print(f"[{datetime.now():%H:%M:%S}]  burst detection on folder  {spike_root}")

    #                 folder_path      save_to    id1   id2      method_set             params_set
    experiment_burst_detection(
        spike_root,   save_root,   "csv",  "Spikes",  cfg.detection_methods,  cfg.parameter_dict
    )

    print("✓ burst detection finished — results in", save_root)


# In[ ]:




