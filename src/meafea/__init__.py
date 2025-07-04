#!/usr/bin/env python
"""
meafea package façade
────────────────────────────────────────────────────────
Provides

    from meafea import MEAFEA

`MEAFEA` wraps:
•  io.discover_paths()              – your folder / Excel discovery helper
•  spike_detection.SpikeConfig
•  spike_detection.run_spike_detection_serial
and keeps DataFrame attributes (`spike_table`, `burst_table`, …)
ready for notebooks.
"""

from __future__ import annotations
from importlib import import_module
from pathlib import Path
from typing import Optional, Union, Sequence
import pandas as pd
from meafea.Templates_vs1 import *

# -------------------------------------------------------------------
# import sub-modules (must exist as files in the same package)
# -------------------------------------------------------------------
io = import_module("meafea.io")
sd = import_module("meafea.spike_detector")
SpikeConfig = sd.SpikeConfig
run_spike_detection_serial = sd.run_spike_detection_serial


# -------------------------------------------------------------------
# public re-exports
# -------------------------------------------------------------------
__all__ = ["MEAFEA", "io", "sd"]


# -------------------------------------------------------------------
# MEAFEA façade class
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# MEAFEA façade class  (put this in src/meafea/__init__.py)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# MEAFEA façade class
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# MEAFEA façade class
# -------------------------------------------------------------------
class MEAFEA:
    """
    Lightweight wrapper around:
        • io.discover_paths()          (your folder discovery helper)
        • spike_detection runner       (serial CSV writer)

    Parameters
    ----------
    root : str | Path
        Experiment top-level folder (contains hdf5/, burst/, Excel, …).
    excel_name : str, optional
        Plate-map Excel file if several exist.
    save_root : str | Path, optional
        Where to write <file>Spikes.csv.  Default: <root>/spike_detection_csv
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        root: str | Path,
        groupcol: str='CellLine',
        *,
        excel_name: str | None = None,
        save_root: str | Path | None = None,
    ):
        from pathlib import Path
        import pandas as pd

        self.root = Path(root).expanduser().resolve()
        self.Name=excel_name
        sel.groupcol=groupcol
        self.RATE_MAP = {
        "MFR" : ("spikedata",  "Timestamp [µs]",  "MFR(Hz)"),
        "BFR" : ("burstdata",   "Timestamp [µs]",  "BFR(Hz)"),
        "NBFR": ("networkburstdata", "Starts",     "NBFR(Hz)"),
    }

        # --- discover folder layout via your io helper ----------------
        self.paths = io.discover_paths(self.root, excel_name=excel_name)

        # helper that works for dict *or* attribute objects
        def _get(obj, key, default):
            if isinstance(obj, dict):
                return obj.get(key, default)
            return getattr(obj, key, default)

        # --- folder to hold spike CSVs -------------------------------
        self.save_root = Path(
            save_root if save_root is not None
            else _get(self.paths, "spike_csv",
                      self.root / "spike_detection_csv")
        ).expanduser().resolve()
        self.save_root.mkdir(parents=True, exist_ok=True)

        # cached DataFrames
        self.spikedata: pd.DataFrame | None = None
        self.burstdata: pd.DataFrame | None = None

        # store helper for later use
        self._get = _get
        
        self.template_df=template_excel(self.root, self.Name)
        self.excelpath, self.excel=find_excel(self.root, self.Name)

    # ----------------------------------------------------------------
    # spike-detection wrapper
    # ----------------------------------------------------------------
    def detect_spikes(self, **cfg_overrides):
        """
        Run spike detection once for *all* hdf5 recordings.

        Keyword arguments override any SpikeConfig field.
        """
        meta_mask = cfg_overrides.pop("MetadataMaskPath", None)

        cfg = SpikeConfig(MetadataMaskPath = meta_mask or str(
            self._get(self.paths, "metadata_mask",
                      self.root / "Metadata")),
        **cfg_overrides,          # now MetadataMaskPath is NOT duplicated
    )
        run_spike_detection_serial(
            hdf5_root=self._get(self.paths, "hdf5", self.root / "hdf5"),
            save_root=self.save_root,
            cfg=cfg,
        )
        

    # ----------------------------------------------------------------
    # convenient loader
    # ----------------------------------------------------------------
    def load_spikes(self, refresh: bool = False):
        """Concatenate every *Spikes.csv* in `self.save_root`."""
        import pandas as pd
        if self.spikedata is None or refresh:
            csvs = list(self.save_root.glob("*Spikes.csv"))
            if not csvs:
                raise FileNotFoundError(
                    "No Spikes.csv found; run detect_spikes() first."
                )
            self.spikedata = pd.concat(pd.read_csv(f) for f in csvs)
        

    # ----------------------------------------------------------------
    # pretty representation
    # ----------------------------------------------------------------
    def __repr__(self):
        n = len(self.spike_table) if self.spike_table is not None else 0
        return f"<MEAFEA root='{self.root.name}'  spikes={n} rows>"
    
    
    def detect_bursts(self, **cfg_kwargs):
        """
        Folder-wide burst detection.

        Any keyword here is passed into BurstConfig, e.g.
            methods=['maxinterval'], parameter_dict={'maxinterval': {...}}
        """
        from meafea.burst_detector import BurstConfig, run_burst_detection_serial
        
        burst_dir = self._get(self.paths, "burst", self.root / "burst")
        burst_dir = Path(burst_dir).expanduser().resolve()
        self.burst_dir=burst_dir
        burst_dir.mkdir(parents=True, exist_ok=True)

        cfg = BurstConfig(
            detection_methods = cfg_kwargs.get("detection_methods",
                                               ["maxinterval"]),
            parameter_dict    = cfg_kwargs.get("parameter_dict",
                                               {"maxinterval": {"set7": [30000, 40000, 100000, 5, 100000]}}),
        )

        run_burst_detection_serial(
            spike_root = self.save_root,            # where Spikes.csv live
            save_root  = burst_dir,
            cfg        = cfg,)

        print("✓ burst detection finished.")
        
    def load_bursts(self, refresh: bool = False):
        import pandas as pd
        if self.burstdata is None or refresh:
            csvs = list(self.burst_dir.glob("*Bursts"))
            if not csvs:
                raise FileNotFoundError(
                    "No Bursts found; run detect_bursts() first."
                )
            self.burstdata = pd.concat(pd.read_csv(f) for f in csvs)
            
            
    def detect_network_burst(self, params_network=2):
        
        #groupby=['Well ID', 'Well Label', 'Dose Label',
       #'Compound ID', 'Experiment','File','Group', 'Cells', 'Culture', 'Label', 'LEDon']
        
        datab=self.burstdata

        network_burst=datab.groupby(self.bywell).apply(lambda x: network(x, params_network)).reset_index()
        
        self.networkburstdata=network_burst
        
    def groupies(self, add_columns: list[str] | None = None) -> None:
    """
    Populate the grouping-attribute shortcuts used throughout MEAFEA.

    Parameters
    ----------
    add_columns : list[str] | None, optional
        Extra column names to append to every base grouping.  Defaults
        to an empty list (safe — no mutable default).
    """
    add_columns = list(dict.fromkeys(add_columns or []))  # deduplicate

    # ---- base groups ------------------------------------------------
    bychannel_base = [
        "Channel ID", "Well ID", "Well Label", "Channel Label",
        "Experiment", "Dose Label", "Compound ID", "File",
        "True Label", "Group",
    ]

    bywell_base = [
        "Well ID", "Well Label", "Experiment", "Dose Label",
        "Compound ID", "File", "Group",
    ]

    bywellnof_base = [
        "Well ID", "Well Label", "Experiment", "Dose Label",
        "Compound ID", "Group",
    ]

    bychannelnof_base = [
        "Channel ID", "Well ID", "Well Label", "Experiment",
        "Dose Label", "Compound ID", "Group",
    ]

    # ---- public attributes -----------------------------------------
    self.bychannel     = bychannel_base     + add_columns
    self.bywell        = bywell_base        + add_columns
    self.bywellnof     = bywellnof_base     + add_columns
    self.bychannelnof  = bychannelnof_base  + add_columns
    self.bygroup=['Experiment', 'Dose Label', 'Compound ID',  'File', 'Group']
    self.bygroupnof=['Experiment', 'Dose Label', 'Compound ID']
  
    def feed_template(self, attrib: str) -> None:
        """
        Merge the object’s *template_df* into the DataFrame stored in
        ``self.<attrib>`` (spike / burst / etc.), adding missing wells or
        channels.  Result overwrites the attribute.

        Works for both “channel mode” (data have *Channel ID*) and “well mode”.
        """
        import pandas as pd

        df = getattr(self, attrib).copy()
        grp_cols = self.bygroupnof          # experiment / dose / compound …

        channel_mode = "Channel ID" in df.columns

        # ----------------------------------------------------------------
        # choose merge keys & template slice
        # ----------------------------------------------------------------
        if channel_mode:
            merge_on = ["Channel ID", "Well Label", "Well ID"]
            template  = self.template_df
            df.drop(columns=["Channel Label"], errors="ignore", inplace=True)
        else:
            merge_on = ["Well ID", "Well Label"]
            template = (self.template_df
                        .drop_duplicates("Well ID")
                        .drop(columns=["Channel ID", "Channel Label"], errors="ignore"))

        # columns that are NOT grouping keys (stay unique per row)
        payload_cols = [c for c in df.columns if c not in grp_cols]

        # ----------------------------------------------------------------
        # single, vectorised merge  (much faster than per-group apply)
        # ----------------------------------------------------------------
        merged = (
            df[grp_cols + payload_cols]
            .merge(template, on=merge_on, how="outer")
            .fillna(0)
        )

        # ensure uniqueness per group + ID
        dedup_key = grp_cols + (["Channel ID"] if channel_mode else ["Well ID"])
        merged = merged.drop_duplicates(subset=dedup_key, keep="first")

        # ----------------------------------------------------------------
        # store back
        # ----------------------------------------------------------------
        setattr(self, attrib, merged)
    

    # ------------------------------------------------------------------
    def load_meta(self, *, layoutmethod: str = "npy", Layout: bool = True) -> None:
        """
        Load stamp / LED / Lab-Book metadata from *self.spikepath* and
        optionally attach the well-layout array.

        Parameters
        ----------
        layoutmethod : {"npy", "excel"}, default "npy"
            Where `find_plate_layout()` should look for the grouping array.
        Layout : bool, default True
            If True and a layout is found, store it in `self.layout`;
            otherwise fallback to 24 × "Cells".
        """
        dpath = Path(self.save_root).expanduser().resolve()

        def _safe_load(name: str):
            f = dpath / f"{name}.npy"
            if f.exists():
                return np.load(f, allow_pickle=True).item()
            print(f"[WARN] {f.name} not found in {dpath}")
            return {}

        # --- core metadata ------------------------------------------------
        self.stamp = _safe_load("stamp")
        self.LED   = _safe_load("LED_info")
        self.Book  = _safe_load("Lab_Book")

        # --- plate layout -------------------------------------------------
        group_arr = self.find_plate_layout(dpath, method=layoutmethod)

        if Layout and len(group_arr):
            self.layout = group_arr
        else:
            self.layout = np.repeat(["Cells"], 24)


    def feed_benchling(self, Param):

            self.template_df=Benchling_layout(self.root, self.template_df, Param)

    def find_plate_layout(self,

                          method: str = "npy"):
        """
        Return the well–group array (`Group_iN`) for an experiment folder.

        Parameters
        ----------
        dpathexperiment : str | Path
            Folder containing either *layout*.npy or the Excel that was
            already parsed into `self.template_df`.
        method : {"npy", "excel"}, default "npy"
            • "npy"   – load first file whose name contains "layout"
            • "excel" – derive groups from `self.template_df`

        Notes
        -----
        Falls back to an empty list and prints a hint if no layout is found.
        """
        dpath = self.root
        group_in: list | np.ndarray = []

        # -------------------------------------------------------------- npy
        if method.lower() == "npy":
            try:
                layout_path = next(p for p in dpath.iterdir() if "layout" in p.name)
                group_in = np.load(layout_path, allow_pickle=True)
            except StopIteration:
                print("[WARN] no *layout*.npy file found in", dpath)
            except Exception as e:
                print("[ERROR] could not load", layout_path.name, "→", e)

        # -----------------------------------------------------------  excel
        elif method.lower() == "excel":
            map_excel = getattr(self, "template_df", None)

            if map_excel is not None and not map_excel.empty:
                merged = map_excel.drop_duplicates("Well ID", keep="last")
                group_in = (
                    merged.sort_values("Well ID")[self.groupcol]  # 'Group' col
                    .values
                )
            else:
                print("[WARN] template_df is empty – create / load Excel first")

        else:
            print("[WARN] unknown method:", method)

        return group_in

    # ------------------------------------------------------------------
    # 1)  Generic threshold filter (writes new attr “…Filtered<param><val>”)
    # ------------------------------------------------------------------
    def filter_by_dict(self, spec: dict[str, dict[str, float]]) -> None:
        """
        spec = {
            "spikedata": {"MFR(Hz)": 0.5},
            "burstdata": {"Spike Count": 30},
            ...
        }
        Creates new attributes called e.g.  spikedataFilteredMFR(Hz)0.5
        """
        for df_attr, cond in spec.items():
            param, value = next(iter(cond.items()))
            df = getattr(self, df_attr)
            mask = df[param] >= value

            # optional grouping: keep group structure but without slow apply
            grouped = df.groupby(self.bygroupnof, group_keys=False)
            filtered = grouped.filter(lambda g: (g[param] >= value).any())

            setattr(self, f"{df_attr}Filtered{param}{value}", filtered)


    # ------------------------------------------------------------------
    # 2)  Filter by list of pre-selected channels
    # ------------------------------------------------------------------
    def filter_by_selected(self, datasets: list[str]) -> None:
        """
        Keep only rows whose 'Channel ID' is in self.selectedchannels.
        Writes <attr>Filtered back to the object.
        """
        for df_attr in datasets:
            df = getattr(self, df_attr)
            mask = df["Channel ID"].isin(self.selectedchannels)
            setattr(self, f"{df_attr}Filtered", df[mask])


    # ------------------------------------------------------------------
    # 3)  Channel-selection helper (creates self.selectedchannels)
    # ------------------------------------------------------------------
    from functools import reduce


    def select_channels_by_dict(
        self,
        spec: dict[str, dict[str, float]],
        doses: list[str] = ("Control",),
    ) -> None:
        """
        spec = {"spikedata": {"MFR(Hz)": 0.5}, "burstdata": {"Spike Count": 30}}
        Result stored in self.selectedchannels (∩ across all rules).
        """
        selections = {}
        for df_attr, cond in spec.items():
            param, value = next(iter(cond.items()))
            df = getattr(self, df_attr)
            df = df[df["Dose Label"].isin(doses)]
            selections[param] = df.loc[df[param] >= value, "Channel ID"].unique()

        self.selectedchannels = reduce(np.intersect1d, selections.values())
        self.selection_dict = selections


    # ------------------------------------------------------------------
    # 4)  Append 'Age' column (vectorised)
    # ------------------------------------------------------------------
    def add_age(self, df_attr: str, unit: str = "DIV") -> None:
        df = getattr(self, df_attr).copy()
        df["Age"] = df["Experiment"].map(lambda x: age(x, unit))
        setattr(self, df_attr, df)


    import pickle, tempfile, os

    def save_object(self, filename: str) -> None:
        """Pickle *self* atomically to <savepath>/<filename>.pkl."""
        target = Path(self.root) / f"{filename}.pkl"
        tmp    = target.with_suffix(".tmp")

        with open(tmp, "wb") as fh:
            pickle.dump(self, fh, pickle.HIGHEST_PROTOCOL)

        os.replace(tmp, target)         # atomic move on the same drive


    # helper (reuse in both methods)


    def Mean_Rates(self, name: str) -> None:
        attr, ts_col, label = self.RATE_MAP[name]
        df = getattr(self, attr)

        out = (
            df.groupby(self.bychannel if name != "NBFR" else self.bywell)
              .apply(lambda x: calcspontanmfr(x, self.LED, self.stamp, ts_col))
              .reset_index(name=f"{label}(Hz)")
        )
        setattr(self, label, out)


    def Mean_Rates_LED_Post(self, name: str, *, post: bool = True) -> None:
        attr, ts_col, label = self.RATE_MAP[name]
        df = getattr(self, attr)

        out = (
            df.groupby(self.bychannel if name != "NBFR" else self.bywell)
              .apply(lambda x: calcspontanmfr(x, self.LED, self.stamp,
                                              ts_col, post=post))
              .reset_index(name=f"{label}(Hz)")
        )
        setattr(self, f"{label}Postled", out)
        
    from functools import reduce

    def StatofParams(self,
                     data_attr: str,
                     params: list[str],
                     by: list[str],
                     *,
                     name: str = "Channel") -> None:
        df = getattr(self, data_attr)

        stats = [
            df.groupby(by)
              .apply(lambda x: param_DistStat(x, p))
              .reset_index()
              .rename(columns=lambda c: c if c in by else f"{p}_{c}")
            for p in params
        ]

        merged = reduce(
            lambda a, b: pd.merge(a, b, how="outer", on=by),
            stats
        ).fillna(0)

        setattr(self, f"{data_attr}statby{name}", merged)

    def CombineParams(self,
                      src_attrs: list[str],
                      by: list[str],
                      new_attr: str) -> None:

        by = [c for c in by if c != "File"]

        dfs = []
        for attr in src_attrs:
            df = getattr(self, attr)
            if "File" in df.columns:
                df = df.drop(columns="File")
            dfs.append(df)

        combined = reduce(lambda a, b: pd.merge(a, b, on=by, how="outer"), dfs)
        setattr(self, new_attr, combined)

# ------------------------------------------------------------------
# 1)  Inter-event intervals (ISI / IBI …)


    def InterIntervals(self,
                       by: list[str],
                       data_attr: str,
                       ts_col: str,
                       name: str) -> None:
        """
        Compute irregular-interval stats on *ts_col* and store in
        self.<data_attr><name>.
        """
        df = getattr(self, data_attr)

        out = (
            df.groupby(by, group_keys=False)
              .apply(lambda g: ISIStat(g, ts_col, name))
              .reset_index()
        )
        setattr(self, f"{data_attr}{name}", out)


    # ------------------------------------------------------------------
    # 2)  Add simple burst / network-burst features
    # ------------------------------------------------------------------
    def More_BurstFeatures(self) -> None:
        bd = self.burstdata
        bd["Burst_Median_Inter-Spike_Interval"] = bd["Spikes"].map(np.median)
        bd["Burst_Median/Mean_ISI"]             = bd["Spikes"].map(
            lambda x: np.median(x) / np.mean(x) if len(x) else 0
        )

    def More_NetworkBurstFeatures(self) -> None:
        nbd = self.networkburstdata
        nbd["Median/Mean_ISI_NetworkBurst"] = nbd["Spikes All"].map(
            lambda x: np.median(x) / np.mean(x) if len(x) else 0
        )
        nbd["Number_elecs_NetworkBurst"] = nbd["Channel Labels"].str.len()
        nbd["Number_spikes_perelec_NetworkBurst_avg"] = nbd["Spikes per channel"].map(
            lambda lst: np.mean([len(sp) for sp in lst]) if lst else 0
        )
        nbd["Number_spikes_perelec_NetworkBurst_std"] = nbd["Spikes per channel"].map(
            lambda lst: np.std([len(sp) for sp in lst]) if lst else 0
        )
        nbd["Start_electrode"] = nbd["Channel Labels"].str[0]


    # ------------------------------------------------------------------
    # 3)  Column-definition helpers (unchanged API, clearer layout)
    # ------------------------------------------------------------------
    def BurstCols(self) -> None:
        self.burscols = [
            "Duration", "Spike Count", "FT Spikes", "Max ISI", "Min ISI",
            "Mean ISI", "Variance", "Burstiness",
        ]

        self.changeburstcols = [
            "Burst Duration", "Number of Spikes per Burst",
            "Proportion of single burst spikes", "MaxISI per Burst",
            "MinISI per Burst", "MeanISI per Burst",
            "Variance of ISI per Burst", "Burstiness",
        ]

        self.column_rename_burstdata = {
            "Duration"   : "Burst_Duration_(sec)",
            "Spike Count": "Number_of_Spikes_per_Burst",
            "FT Spikes"  : "First_to_Threshold_Spikes",
            "Max ISI"    : "Burst_Max_Inter-Spike_Interval",
            "Min ISI"    : "Burst_Min_Inter-Spike_Interval",
            "Mean ISI"   : "Burst_Mean_Inter-Spike_Interval",
            "Variance"   : "Burst_Variance_of_Spike_Times",
            "Burstiness" : "Burstiness_Index",
        }
        self.newburstcols = list(self.column_rename_burstdata.values())


    def NetworkBurstCols(self) -> None:
        self.allnetworkburstcols = [
            "SpikeLength", "MinISI", "MaxISI", "MeanxISI", "MedianxISI",
            "Peakslen", "PeaksValmax", "PeakProminence", "PeakWidth", "MaxProm",
            "MaxPeakWidth", "NetworkMeanSTTC", "NetworkMedianSTTC",
            "Network25STTC",
        ]

        self.modnetworkburstcols = [
            "Number of Spikes per Network Burst", "MinISI per NB",
            "MaxISI per NB", "MeanxISI per NB", "MedianISI per NB",
            "N of peak in FR per NB", "Max FR per NB", "Slope of FR peak per NB",
            "Peak Duration per NB", "Max Slope per NB", "Max FR peak per NB",
            "Mean STTC per NB", "Median STTC per NB", "Q1 STTC per NB",
        ]

        self.column_rename_networkburstdata = {
            "SpikeLength"       : "Network_Spike_Length_(ms)",
            "MinISI"            : "Network_Min_Inter-Spike_Interval",
            "MaxISI"            : "Network_Max_Inter-Spike_Interval",
            "MeanxISI"          : "Network-Mean_Inter-Spike_Interval",
            "MedianxISI"        : "Network_Median_Inter-Spike_Interval",
            "Peakslen"          : "Network_Number_of_Peaks",
            "PeaksValmax"       : "Network_Max_Peak_Value",
            "PeakProminence"    : "Network_Peak_Prominence",
            "PeakWidth"         : "Network_Peak_Width",
            "MaxProm"           : "Network_Max_Peak_Prominence",
            "MaxPeakWidth"      : "Network_Max_Peak_Width",
            "NetworkMeanSTTC"   : "Network_Mean_Spike_Time_Tiling_Coefficient",
            "NetworkMedianSTTC" : "Network_Median_Spike_Time_Tiling_Coefficient",
            "Network25STTC"     : "Network_25th_Percentile_Spike_Time_Tiling_Coefficient",
        }
        self.newnetworkburstcols = list(self.column_rename_networkburstdata.values())


    # ------------------------------------------------------------------
    # 4)  Generic renamer (minor fix)
    # ------------------------------------------------------------------
    def rename_columns(self, data_attr: str) -> None:
        mapping = getattr(self, f"column_rename_{data_attr}")
        df = getattr(self, data_attr).rename(columns=mapping)
        setattr(self, data_attr, df)


    # ------------------------------------------------------------------
    # 5)  Lists of non-parameter columns
    # ------------------------------------------------------------------
    def NonParamColumns(self) -> None:
        self.MEANonParamcolumnsByWell = [
            "Experiment", "Dose Label", "Compound ID", "Well ID",
            "Well Label", "Group", "Age",
        ]

        self.MEANonParamcolumnsByChannel = [
            "Experiment", "Dose Label", "Compound ID", "Well ID", "Well Label",
            "Channel ID", "True Label", "Channel Label", "Group", "Age",
        ]

        self.ExcelNonParamcolumns = [
            "Compound", "Concentration_µM", "Vehicle", "Treatment",
            "Experiment_Excel", "hdf5_data_folder", "ra# ------------------------------------------------------------------
#  GUI helpers – lazy import + simple guard
# ------------------------------------------------------------------
    def ViewRaw(self) -> None:
        """
        Launch the Tkinter viewer (`MEATkinterApp`) for raw files.
        """
        try:
            from meafea.MEATkinter_class_vs1 import MEATkinterApp
            MEATkinterApp(self, "Dose Label").tkstart()
        except Exception as exc:
            print("[ViewRaw] GUI could not start →", exc)


    def ChooseThreshold(self) -> None:
        """
        Launch the interactive threshold-selection GUI (`MEAThreshold`).
        """
        try:
            from meafea.MEATkinter_class_vs1 import MEAThreshold
            MEAThreshold(self).tkstart()
        except Exception as exc:
            print("[ChooseThreshold] GUI could not start →", exc)
    w_data_folder",
                "save_folder", "PlateID", "Date_neuronal_plating",
                "Date_compound_treatment", "Cell_origin", "Genotype", "CellLine",
                "Cell_culture_type", "Neuron_density_k_per_well",
                "Astrocyte_density_k_per_well", "Media_type", "Comment", "Well",
            ]

