"""
Ultrasonic NDT pulse-echo acquisition (Functional Approach).
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Iterable, List, Optional, Dict, Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import uniform_filter1d
# --------------------------------------------------------------------------- #
#  Device handle                                                              #
# --------------------------------------------------------------------------- #

_PROBE = None


def get_probe(force_new: bool = False):
    global _PROBE
    if _PROBE is None or force_new:
        from lib.device import Pic0rick
        _PROBE = Pic0rick()
    return _PROBE


# --------------------------------------------------------------------------- #
#  Core Data Helpers                                                          #
# --------------------------------------------------------------------------- #

def get_time_vector(acq: Dict[str, Any]) -> np.ndarray:
    return np.arange(len(acq["signal"])) / acq["Fech"]


def get_key(acq: Dict[str, Any]) -> tuple:
    return (
        round(float(acq["gain"]), 6),
        int(acq["pon"]),
        int(acq["poff"]),
        int(acq["damp"]),
        str(acq.get("target", "")),
        str(acq.get("piezo_id", "")),
    )


def get_label(acq: Dict[str, Any]) -> str:
    target = acq.get("target") or "untitled"
    parts = [
        target, 
        f"gain={acq['gain']} dB",
        f"pon={acq['pon']}, poff={acq['poff']}, damp={acq['damp']}",
        f"Fech={acq['Fech'] / 1e6:.2f} MHz"
    ]
    if acq.get("piezo_id"):
        parts.insert(1, f"piezo={acq['piezo_id']}")
    return "  |  ".join(parts)


def filter_signal(acq: Dict[str, Any]) -> np.ndarray:
    freq = acq.get("piezo_central_freq")
    bw = acq.get("piezo_bandwidth")
    if freq is None or bw is None:
        raise ValueError("Filtering needs piezo_central_freq and piezo_bandwidth.")
    nyq = acq["Fech"] / 2.0
    low = max((freq - bw / 2) / nyq, 1e-6)
    high = min((freq + bw / 2) / nyq, 1 - 1e-6)
    sos = butter(4, [low, high], btype="bandpass", output="sos")
    return sosfiltfilt(sos, acq["signal"]).astype(np.float32)


def get_amplitude(acq: Dict[str, Any], start_us: float, end_us: float, filtered: bool = True) -> float:
    sig = filter_signal(acq) if filtered else acq["signal"]
    t_us = get_time_vector(acq) * 1e6
    mask = (t_us >= start_us) & (t_us <= end_us)
    window = sig[mask]
    return float(np.max(np.abs(window))) if window.size else 0.0


# --------------------------------------------------------------------------- #
#  Acquisition & Processing                                                   #
# --------------------------------------------------------------------------- #

def acquire_signal(
    Fech: float,
    *,
    probe=None,
    pon: int = 70,
    poff: int = 70,
    damp: int = 6000,
    gain: float = 20,
    target: str = "",
    piezo_id: str = "",
    piezo_central_freq: Optional[float] = None,
    piezo_bandwidth: Optional[float] = None,
    h5_path: Optional[str] = None,
    overwrite: bool = False,
) -> Dict[str, Any]:
    
    if h5_path is not None and not overwrite:
        key = (round(float(gain), 6), int(pon), int(poff), int(damp), str(target), str(piezo_id))
        cached = _load_by_key(h5_path, key)
        if cached is not None:
            dirty = False
            if piezo_central_freq is not None and cached.get("piezo_central_freq") != piezo_central_freq:
                cached["piezo_central_freq"] = piezo_central_freq
                dirty = True
            if piezo_bandwidth is not None and cached.get("piezo_bandwidth") != piezo_bandwidth:
                cached["piezo_bandwidth"] = piezo_bandwidth
                dirty = True
            if dirty:
                save_acquisition(cached, h5_path)
            return cached

    if probe is None:
        probe = get_probe()
    probe.dac(gain)
    probe.pulse_adc_trigger(pon=pon, poff=poff, damp=damp)
    C = probe.read()
    raw = [x.replace("b'", "") for x in str(C[2]).split(",") if len(x)]
    signal = np.array([(int(x, 16) - 512) / 512.0 for x in raw[:-1]], dtype=np.float32)
    
    acq = {
        "signal": signal,
        "Fech": Fech,
        "pon": pon,
        "poff": poff,
        "damp": damp,
        "gain": gain,
        "target": target,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "piezo_id": piezo_id,
        "piezo_central_freq": piezo_central_freq,
        "piezo_bandwidth": piezo_bandwidth,
        "h5_path": h5_path
    }
    
    if h5_path is not None:
        save_acquisition(acq, h5_path)
    return acq


def detect_echoes(
    acq: Dict[str, Any],
    *,
    start_us: float,
    end_us: float,
    target_thickness: float,
    speed_of_sound: float,
    smooth: int = 10,
    smooth_kernel: int = 5,
    tolerance: float = 0.30,
    min_amplitude_ratio: float = 0.10,
    dip_threshold: float = 0.5,
    plot: bool = True,
) -> dict:
    filtered = filter_signal(acq)
    work = np.abs(filtered)
    for _ in range(int(smooth)):
        work = uniform_filter1d(work, size=int(smooth_kernel), mode="nearest")

    t_us = get_time_vector(acq) * 1e6
    dt_us = 2.0 * target_thickness / speed_of_sound * 1e6

    main_mask = (t_us >= start_us) & (t_us <= end_us)
    if not np.any(main_mask):
        raise ValueError(f"Window [{start_us}, {end_us}] µs is outside the signal (0..{t_us[-1]:.2f} µs).")
    
    main_win = np.where(main_mask)[0]
    main_idx = int(main_win[np.argmax(work[main_win])])
    main_amp = float(work[main_idx])
    main_time_us = float(t_us[main_idx])

    peak_indices = [main_idx]
    peak_amps = [main_amp]
    peak_n = [0]
    thresh = min_amplitude_ratio * main_amp
    tol_us = tolerance * dt_us
    n = 1
    
    while True:
        center = main_time_us + n * dt_us
        if center > end_us:
            break
        lo = max(center - tol_us, start_us)
        hi = min(center + tol_us, end_us)
        win = np.where((t_us >= lo) & (t_us <= hi))[0]
        if win.size == 0:
            n += 1
            continue
        cand_idx = int(win[np.argmax(work[win])])
        cand_amp = float(work[cand_idx])
        if cand_amp >= thresh:
            peak_indices.append(cand_idx)
            peak_amps.append(cand_amp)
            peak_n.append(n)
        n += 1

    peak_times_us = [float(t_us[i]) for i in peak_indices]
    discarded_times_us: List[float] = []
    
    if len(peak_amps) >= 3:
        keep = [True] * len(peak_amps)
        for i in range(1, len(peak_amps) - 1):
            neighbour_min = min(peak_amps[i - 1], peak_amps[i + 1])
            if peak_amps[i] < dip_threshold * neighbour_min:
                keep[i] = False
        discarded_times_us = [peak_times_us[i] for i, k in enumerate(keep) if not k]
        peak_indices  = [peak_indices[i] for i, k in enumerate(keep) if k]
        peak_amps     = [peak_amps[i] for i, k in enumerate(keep) if k]
        peak_n        = [peak_n[i] for i, k in enumerate(keep) if k]
        peak_times_us = [peak_times_us[i] for i, k in enumerate(keep) if k]

    if len(peak_times_us) >= 2:
        first_us = peak_times_us[0]
        last_us  = peak_times_us[-1]
        n_intervals = peak_n[-1] - peak_n[0]
        measured_dt_us = (last_us - first_us) / n_intervals
        calculated_thickness = measured_dt_us * 1e-6 * speed_of_sound / 2.0
    else:
        measured_dt_us = None
        calculated_thickness = None
        n_intervals = 0

    fig = None
    if plot:
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(t_us, acq["signal"], lw=0.8, color="tab:blue", alpha=0.45, label="raw")
        ax.plot(t_us, work, lw=1.1, color="tab:orange", label="|filtered| smoothed")
        ax.axvspan(start_us, end_us, color="gray", alpha=0.10, label=f"window [{start_us}, {end_us}] µs")
        for k, tp in enumerate(peak_times_us):
            ax.axvline(tp, color="tab:red", lw=1.0, ls="-" if k == 0 else "--", alpha=0.85)
        for tp in discarded_times_us:
            ax.axvline(tp, color="tab:gray", lw=1.0, ls=":", alpha=0.6)

        ax.plot([], [], color="tab:red", lw=1.0, ls="--", label=f"{len(peak_times_us)} peak(s) kept")
        if discarded_times_us:
            ax.plot([], [], color="tab:gray", lw=1.0, ls=":", label=f"{len(discarded_times_us)} dip(s) discarded")

        title_lines = [
            get_label(acq),
            f"expected Δt = {dt_us:.2f} µs (thickness={target_thickness*1e3:.2f} mm, c={speed_of_sound:.0f} m/s)",
        ]
        if calculated_thickness is not None:
            title_lines.append(f"measured Δt = {measured_dt_us:.2f} µs → calculated thickness = {calculated_thickness*1e3:.3f} mm")
        ax.set_title("\n".join(title_lines))
        ax.set_xlabel("Time [µs]")
        ax.set_ylabel("Amplitude [norm.]")
        ax.grid(True, alpha=0.3)
        ax.margins(x=0)
        ax.legend(loc="upper right", fontsize=9)

    return {
        "filtered_signal": filtered,
        "work_signal": work,
        "peak_indices": peak_indices,
        "peak_times_us": peak_times_us,
        "peak_amplitudes": peak_amps,
        "discarded_times_us": discarded_times_us,
        "main_peak_time_us": main_time_us,
        "expected_dt_us": dt_us,
        "measured_dt_us": measured_dt_us,
        "calculated_thickness": calculated_thickness,
        "figure": fig,
    }


def calibrate(
    h5_path: str,
    *,
    Fech: float,
    start_us: float,
    end_us: float,
    gain: float,
    piezo_central_freq: float,
    piezo_bandwidth: float,
    piezo_id: str = "",
    target: str = "calibration",
    damp: int = 6000,
    pon_poff_values: Iterable[int] = range(25, 156, 10),
    probe=None,
    overwrite: bool = False,
) -> dict:
    values = list(pon_poff_values)
    amps = np.empty(len(values))
    acqs = []

    for i, p in enumerate(values):
        acq = acquire_signal(
            Fech=Fech, probe=probe, gain=gain, pon=p, poff=p, damp=damp, 
            target=target, piezo_id=piezo_id, piezo_central_freq=piezo_central_freq,
            piezo_bandwidth=piezo_bandwidth, h5_path=h5_path, overwrite=overwrite,
        )
        amps[i] = get_amplitude(acq, start_us, end_us)
        acqs.append(acq)

    i_max, i_min = int(np.argmax(amps)), int(np.argmin(amps))

    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2, 1], hspace=0.35, wspace=0.25)
    ax_min, ax_max = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax_amp = fig.add_subplot(gs[:, 1])

    for ax, idx, tag in ((ax_min, i_min, "MIN"), (ax_max, i_max, "MAX")):
        acq = acqs[idx]
        t_us = get_time_vector(acq) * 1e6
        ax.plot(t_us, filter_signal(acq), lw=0.8)
        ax.axvspan(start_us, end_us, color="gray", alpha=0.15, label=f"window [{start_us}, {end_us}] µs")
        ax.set_title(f"{tag} — pon = poff = {values[idx]} (amplitude = {amps[idx]:.4f})")
        ax.set_ylabel("Amplitude [norm.]")
        ax.grid(True, alpha=0.3)
        ax.margins(x=0)
    ax_max.set_xlabel("Time [µs]")

    ax_amp.plot(values, amps, "o-", lw=1.2)
    ax_amp.plot(values[i_max], amps[i_max], "o", color="tab:red", markersize=10, label=f"best: pon=poff={values[i_max]}")
    ax_amp.plot(values[i_min], amps[i_min], "o", color="tab:gray", markersize=8, label=f"worst: pon=poff={values[i_min]}")
    ax_amp.set_xlabel("pon = poff")
    ax_amp.set_ylabel("Peak amplitude (filtered) [norm.]")
    ax_amp.set_title("Amplitude vs pulse width")
    ax_amp.grid(True, alpha=0.3)
    ax_amp.legend(loc="best", fontsize=9)
    fig.suptitle(f"Calibration — target='{target}', gain={gain} dB, piezo {piezo_central_freq/1e6:.2f} MHz ± {piezo_bandwidth/2/1e6:.2f} MHz")

    return {
        "best_pon_poff": values[i_max],
        "best_amplitude": float(amps[i_max]),
        "worst_pon_poff": values[i_min],
        "worst_amplitude": float(amps[i_min]),
        "pon_poff_values": values,
        "amplitudes": amps.tolist(),
        "figure": fig,
    }


def plot_acquisition(
    acq: Dict[str, Any],
    ax: Optional[plt.Axes] = None,
    unit: str = "us",
    figsize=(15, 5),
    **plot_kwargs,
) -> plt.Axes:
    has_filter = acq.get("piezo_central_freq") is not None and acq.get("piezo_bandwidth") is not None
    scale = {"s": 1.0, "ms": 1e3, "us": 1e6}[unit]
    t = get_time_vector(acq) * scale

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    lw = plot_kwargs.pop("lw", 0.8)
    ax.plot(t, acq["signal"], lw=lw, color="tab:blue", alpha=0.45, label="raw", **plot_kwargs)
    if has_filter:
        ax.plot(t, filter_signal(acq), lw=lw, color="tab:orange", label="filtered")
        ax.legend(loc="upper right", fontsize=9)

    ax.set_xlabel(f"Time [{unit}]")
    ax.set_ylabel("Amplitude [norm.]")
    ax.set_title(get_label(acq))
    ax.grid(True, alpha=0.3)
    ax.margins(x=0)
    return ax


# --------------------------------------------------------------------------- #
#  HDF5 I/O                                                                   #
# --------------------------------------------------------------------------- #

_SCALAR_ATTRS = ("Fech", "pon", "poff", "damp", "gain", "target", "timestamp", "piezo_id", "piezo_central_freq", "piezo_bandwidth")


def _attr(g, name):
    v = g.attrs[name]
    return v.decode() if isinstance(v, bytes) else v


def _attr_or(g, name, default=None):
    return _attr(g, name) if name in g.attrs else default


def _opt_float(v) -> Optional[float]:
    if v is None: return None
    f = float(v)
    return None if np.isnan(f) else f


def _find_group_by_key(h5, key) -> Optional[str]:
    for name in h5:
        g = h5[name]
        try:
            existing = (
                round(float(_attr(g, "gain")), 6),
                int(_attr(g, "pon")),
                int(_attr(g, "poff")),
                int(_attr(g, "damp")),
                str(_attr(g, "target")),
                str(_attr_or(g, "piezo_id", "")),
            )
        except KeyError:
            continue
        if existing == key:
            return name
    return None


def save_acquisition(acq: Dict[str, Any], path: Optional[str] = None) -> None:
    p = path or acq.get("h5_path")
    if p is None:
        raise ValueError("No h5_path provided.")
    acq["h5_path"] = p
    with h5py.File(p, "a") as h5:
        match = _find_group_by_key(h5, get_key(acq))
        if match is not None:
            del h5[match]
        i = 0
        while f"acq_{i:04d}" in h5:
            i += 1
        g = h5.create_group(f"acq_{i:04d}")
        g.create_dataset("signal", data=np.asarray(acq["signal"]), compression="gzip")
        for k in _SCALAR_ATTRS:
            v = acq.get(k)
            if v is not None:
                g.attrs[k] = v


def _load_by_key(path: str, key: tuple) -> Optional[Dict[str, Any]]:
    try:
        h5 = h5py.File(path, "r")
    except (FileNotFoundError, OSError):
        return None
    with h5:
        name = _find_group_by_key(h5, key)
        if name is None: return None
        g = h5[name]
        return {
            "signal": g["signal"][:],
            "Fech": float(_attr(g, "Fech")),
            "pon": int(_attr(g, "pon")),
            "poff": int(_attr(g, "poff")),
            "damp": int(_attr(g, "damp")),
            "gain": float(_attr(g, "gain")),
            "target": str(_attr(g, "target")),
            "timestamp": str(_attr(g, "timestamp")),
            "h5_path": path,
            "piezo_id": str(_attr_or(g, "piezo_id", "")),
            "piezo_central_freq": _opt_float(_attr_or(g, "piezo_central_freq")),
            "piezo_bandwidth": _opt_float(_attr_or(g, "piezo_bandwidth"))
        }


def get_info(h5_path: str, verbose: bool = True) -> dict:
    try:
        h5 = h5py.File(h5_path, "r")
    except (FileNotFoundError, OSError):
        summary = {"h5_path": h5_path, "n_signals": 0, "gains": [], "targets": []}
        if verbose: print(f"HDF5 file not found: {h5_path}")
        return summary
    with h5:
        names = list(h5)
        gains = sorted({float(_attr(h5[n], "gain")) for n in names})
        targets = sorted({str(_attr(h5[n], "target")) for n in names})
    summary = {"h5_path": h5_path, "n_signals": len(names), "gains": gains, "targets": targets}
    if verbose:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


def save_h5(path: str, acquisitions: Iterable[Dict[str, Any]], overwrite_file: bool = False) -> None:
    if overwrite_file:
        with h5py.File(path, "w"): pass
    for acq in acquisitions:
        save_acquisition(acq, path)


def load_h5(path: str) -> List[Dict[str, Any]]:
    out = []
    with h5py.File(path, "r") as h5:
        for name in sorted(h5):
            g = h5[name]
            out.append({
                "signal": g["signal"][:],
                "Fech": float(_attr(g, "Fech")),
                "pon": int(_attr(g, "pon")),
                "poff": int(_attr(g, "poff")),
                "damp": int(_attr(g, "damp")),
                "gain": float(_attr(g, "gain")),
                "target": str(_attr(g, "target")),
                "timestamp": str(_attr(g, "timestamp")),
                "h5_path": path,
                "piezo_id": str(_attr_or(g, "piezo_id", "")),
                "piezo_central_freq": _opt_float(_attr_or(g, "piezo_central_freq")),
                "piezo_bandwidth": _opt_float(_attr_or(g, "piezo_bandwidth"))
            })
    return out
