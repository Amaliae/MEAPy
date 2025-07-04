#!/usr/bin/env python
# coding: utf-8

# In[8]:


from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import pandas as pd
from importlib import import_module


@dataclass
class DiscoveredPaths:
    hdf5: Path
    raw: Path
    save: Path
    spike_dir: Path
    burst_dir: Path
    excel_meta: Optional[pd.DataFrame]
    template_df: Optional[pd.DataFrame]

def discover_paths(
    root: str | Path,
    excel_name: Optional[str] = None,
    identifier: Optional[str] = None,
    nonidentifier: Optional[str] = None,
) -> DiscoveredPaths:
    root = Path(root).expanduser().resolve()
    hdf5, raw, save = root / "hdf5", root, root

    excel_df = template_df = None

    # ---------------- Excel overrides --------------------
    if excel_name:
        excel_path = root / excel_name
        if excel_path.is_file():
            # let pandas figure out which row has the header
            excel_df = pd.read_excel(excel_path, engine="openpyxl", header=0)
            template_df = excel_df.copy()

            def first_value(df: pd.DataFrame, col: str) -> Optional[str]:
                """Return first non-blank value in df[col] or None."""
                if col in df.columns:
                    ser = df[col].dropna().astype(str).str.strip()
                    return ser.iloc[0] if not ser.empty else None
                return None

            for col, current in [
                ("hdf5_data_folder", hdf5),
                ("raw_data_folder", raw),
                ("save_folder", save),
            ]:
                val = first_value(excel_df, col)
                if val:
                    current_path = Path(val).expanduser().resolve()
                    if col == "hdf5_data_folder":
                        hdf5 = current_path
                    elif col == "raw_data_folder":
                        raw = current_path
                    else:
                        save = current_path

    # ---------------- Fallbacks ---------------------------
    if not hdf5.is_dir():
        for p in root.rglob("hdf5"):
            if p.is_dir():
                hdf5 = p
                break
        else:
            raise FileNotFoundError("No 'hdf5' directory under root")

    def _find(folder_name: str) -> Path:
        for p in hdf5.parent.rglob(folder_name):
            if identifier and identifier not in p.name:
                continue
            if nonidentifier and nonidentifier in p.name:
                continue
            return p
        return hdf5.parent / folder_name  # fallback if nothing found

    spike_dir = _find("spike_detection_csv")
    burst_dir = _find("burst")

    (save / "Results").mkdir(parents=True, exist_ok=True)

    return DiscoveredPaths(
        hdf5=hdf5,
        raw=raw,
        save=save,
        spike_dir=spike_dir,
        burst_dir=burst_dir,
        excel_meta=excel_df,
        template_df=template_df,)
 


# In[ ]:





# In[9]:


paths = discover_paths(
    r"E:\AH\NGN2AD2Ratios\iN9\Pharm",
    excel_name="iN9PharmPlatemapMEA.xlsx"    # or None
)

print("HDF5 :", paths.hdf5)
print("Spikes:", paths.spike_dir)
print("Bursts:", paths.burst_dir)


# In[ ]:




