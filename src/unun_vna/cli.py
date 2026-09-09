"""Command-line interface and S-A-A-2 acquisition/calibration code.

Overview
========
The program deliberately separates three jobs:

1. acquire *raw* complex S11/S21 data from the S-A-A-2 via ``pynanovna``;
2. apply a host-side T/R calibration suitable for the S-A-A-2 architecture;
3. hand calibrated data to :mod:`unun_vna.analysis`, which contains the
   fixture and UnUn efficiency mathematics.

Why host-side calibration?
--------------------------
The S-A-A-2 exposes raw receiver data over USB.  ``pynanovna`` can operate the
instrument, but its generic calibration model is not the one we want to rely on
for this measurement.  This tool therefore stores and applies its own
calibration coefficients.

Important architectural limitation
----------------------------------
The S-A-A-2 is a T/R VNA: we measure S11 and forward S21, but not the complete
four-S-parameter set of a bidirectional 2-port VNA.  The host calibration mirrors
the NanoVNA V2 enhanced-response user-calibration path: SOL reflection correction,
SHORT/OPEN leakage removal, THRU normalization using the THRU's own reflection
state, and Port-1 source-match loop-gain correction.  It still cannot perform a
full independent Port-2 load-match correction.

That limitation matters because the one-transformer method often presents a
very high impedance when viewed from Port 2.  The recommended measurement setup
therefore includes a fixed attenuator (typically 10 dB) immediately in front of
Port 2.  The attenuator must remain in place during calibration, fixture
characterization, and DUT measurement.  Its forward loss is calibrated out,
while reflections travelling to/from Port 2 are attenuated twice.

No hidden de-embedding is performed.  All model assumptions are explicit and
the resistor fixture is independently characterized before its data are used
for the efficiency calculation.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import sys
from pathlib import Path

import numpy as np
import pynanovna

from .analysis import (
    ideal_transformer_reference,
    series_fixture_impedance,
    unun_efficiency,
    vswr_from_s11,
)


# Default sweep chosen to cover the HF amateur bands of interest while leaving
# some margin below 40 m and above 10 m.
DEFAULT_START = 1e6
DEFAULT_STOP = 35e6
DEFAULT_POINTS = 501

# Both S-A-A-2 ports use a nominal 50-ohm reference impedance.
DEFAULT_Z0 = 50.0

# Default high-side/low-side impedance ratio. 49 corresponds to an UnUn
# commonly described as 1:49. The ratio is used only for matching/reference
# diagnostics; measured efficiency is reconstructed from measured S-parameters
# and measured fixture impedance.
DEFAULT_RATIO = 49.0

# Frequencies printed by the fixture/analyze commands unless the user overrides
# them with --bands.
DEFAULT_BANDS = "7.1M,14.2M,21.2M,28.5M"


def db20(z):
    """Return the magnitude of a complex S-parameter in dB.

    S-parameters are complex *wave-amplitude* ratios, so magnitude in dB uses
    ``20*log10(|S|)``.  Power ratios later in the efficiency calculation use
    ``|S|**2`` instead.
    """
    a = abs(z)
    return -math.inf if a == 0 else 20.0 * math.log10(a)


def parse_freq(s):
    """Parse a CLI frequency such as ``7100000``, ``7.1M`` or ``7.1MHz``."""
    m = re.fullmatch(
        r"(?i)\s*([0-9]*\.?[0-9]+)\s*([kmg]?)\s*(?:hz)?\s*", s
    )
    if not m:
        raise argparse.ArgumentTypeError(f"bad frequency: {s}")

    # Multipliers are interpreted in the usual SI sense.
    return float(m.group(1)) * {
        "": 1.0,
        "k": 1e3,
        "m": 1e6,
        "g": 1e9,
    }[m.group(2).lower()]


def parse_freq_list(text):
    """Parse a comma-separated list used by ``--bands``."""
    return [parse_freq(x) for x in text.split(",") if x.strip()]


def parse_impedance_ratio(text):
    """Parse the transformer high-side/low-side impedance ratio.

    Accepted examples:

        --ratio 49
        --ratio 1:49
        --ratio 1/49

    All three forms above mean ``Z_high/Z_low = 49``.

    The one-transformer method is normally connected with the lower-impedance
    side at Port 1 and the higher-impedance side toward the fixture/Port 2.
    Ratios below 1 are therefore rejected.
    """
    s = str(text).strip()

    if ":" in s or "/" in s:
        sep = ":" if ":" in s else "/"
        parts = s.split(sep)
        if len(parts) != 2:
            raise argparse.ArgumentTypeError(
                f"bad impedance ratio '{text}'; use e.g. 49 or 1:49"
            )
        try:
            low = float(parts[0])
            high = float(parts[1])
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"bad impedance ratio '{text}'; use e.g. 49 or 1:49"
            ) from exc
        if low <= 0.0 or high <= 0.0:
            raise argparse.ArgumentTypeError("ratio components must be > 0")
        ratio = high / low
    else:
        try:
            ratio = float(s)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"bad impedance ratio '{text}'; use e.g. 49 or 1:49"
            ) from exc

    if not math.isfinite(ratio) or ratio < 1.0:
        raise argparse.ArgumentTypeError(
            "impedance ratio must be >= 1 (high/low); "
            "connect the low-impedance side to Port 1"
        )
    return ratio


def open_vna(index):
    """Open one VNA through pynanovna and verify that it is connected.

    ``pynanovna`` performs the USB/serial discovery.  We keep all interaction
    with that library behind this small helper so that the rest of the program
    deals only with an already-open device object.
    """
    v = pynanovna.VNA(vna_index=index, logging_level="critical")
    if not getattr(v, "connected", False):
        raise RuntimeError("No NanoVNA found")
    return v


def sweep_raw(v, start, stop, points):
    """Acquire one uncalibrated complex S11/S21 sweep.

    The returned arrays are *raw with respect to this program*.  They are not
    yet corrected by ``apply_cal``.

    ``pynanovna`` emits CRITICAL messages saying that no pynanovna calibration
    is active.  That is expected, because this program intentionally uses its
    own host-side calibration.  The messages are suppressed only during this
    raw acquisition call.
    """
    v.set_sweep(start, stop, points)

    previous_disable = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    try:
        s11, s21, f = v.sweep()
    finally:
        logging.disable(previous_disable)

    # Converting once here gives the rest of the program predictable numpy
    # dtypes and allows vectorized complex math.
    return (
        np.asarray(f, float),
        np.asarray(s11, complex),
        np.asarray(s21, complex),
    )


def write_s_csv(path, f, s11, s21):
    """Write complex S11/S21 and their dB magnitudes to a human-readable CSV."""
    with Path(path).open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "frequency_hz",
                "s11_real",
                "s11_imag",
                "s11_db",
                "s21_real",
                "s21_imag",
                "s21_db",
            ]
        )
        for ff, a, b in zip(f, s11, s21):
            w.writerow(
                [
                    int(round(ff)),
                    a.real,
                    a.imag,
                    db20(a),
                    b.real,
                    b.imag,
                    db20(b),
                ]
            )


def capture(v, title, text, start, stop, points):
    """Prompt for a calibration standard and acquire one raw sweep."""
    print("\n" + "=" * 68)
    print(title)
    print(text)
    input("Press Enter when ready...")
    f, s11, s21 = sweep_raw(v, start, stop, points)
    print("Captured.")
    return f, s11, s21


def _apply_reflection_cal(raw11, e00, e11, etrack):
    """Apply the three-term SOL reflection calibration to raw S11.

    The host receives the S-A-A-2 reflection channel after its internal e-cal,
    but before user SOL calibration.  This is the same quantity used by the
    NanoVNA V2 on-screen user-calibration code.
    """
    u = raw11 - e00
    den = etrack + e11 * u
    if np.any(np.abs(den) < 1e-15):
        raise RuntimeError("Singular reflection calibration")
    return u / den


def _source_match_gain(s11, e11):
    """Return the NanoVNA V2 enhanced-response source-match loop gain.

    The firmware calls this quantity ``SOL_compute_thru_gain`` and computes

        gain = 1 / (1 - e11*S11)

    where S11 is already SOL-calibrated.  It describes the signal-flow loop
    between the source-match error of Port 1 and a reflecting DUT.
    """
    den = 1.0 - e11 * s11
    if np.any(np.abs(den) < 1e-15):
        raise RuntimeError("Singular enhanced-response source-match correction")
    return 1.0 / den


def compute_cal(short11, short21, open11, open21, load11, thru11, thru21):
    """Derive the S-A-A-2/NanoVNA V2 enhanced-response T/R calibration.

    This follows the on-screen NanoVNA V2 user-calibration algorithm rather
    than the simpler T/R normalization previously used by this tool.

    Reflection channel
    ------------------
    S11 uses the standard three-term one-port error model::

        m = e00 + etrack*Gamma / (1 - e11*Gamma)

    where ``e00`` is directivity, ``etrack`` is reflection tracking and
    ``e11`` is source match.  SHORT, OPEN and LOAD solve these terms at every
    frequency.

    Transmission leakage
    --------------------
    NanoVNA V2 models forward leakage as an affine function of the *raw*
    reflection channel::

        leakage(raw_s11) = leak0 + leak_r*raw_s11

    The firmware derives this line from the SHORT and OPEN measurements.  Both
    standards ideally have zero true S21, so their measured forward signal is
    treated as leakage/feed-through.

    THRU reference and enhanced response
    ------------------------------------
    The THRU capture must retain both raw S11 and raw S21.  During application
    the THRU transmission is corrected using the THRU's own raw S11, and the
    THRU reference is then scaled by the ratio of source-match loop gains for
    the DUT and THRU.  This is crucial for highly reflecting DUTs such as the
    high-value series fixtures used by the one-transformer method.

    Port-2 load match remains outside this T/R error model.  A fixed attenuator
    at Port 2 is therefore still recommended and must remain in place for both
    calibration and measurement.
    """
    # --- Reflection SOL calibration -----------------------------------------
    e00 = load11
    A = open11 - e00
    B = short11 - e00
    d = A - B
    if np.any(np.abs(d) < 1e-15):
        raise RuntimeError("Degenerate SOL calibration data")
    e11 = (A + B) / d
    etrack = A * (1.0 - e11)

    # --- Forward leakage model ----------------------------------------------
    # Match NanoVNA V2 firmware: use SHORT and OPEN, not LOAD and OPEN.
    d2 = short11 - open11
    if np.any(np.abs(d2) < 1e-15):
        raise RuntimeError("Cannot derive T/R leakage model")
    leak_r = (short21 - open21) / d2
    leak0 = open21 - leak_r * open11

    return e00, e11, etrack, leak0, leak_r, thru11, thru21


def apply_cal(
    raw11,
    raw21,
    e00,
    e11,
    etrack,
    leak0,
    leak_r,
    thru11,
    thru21,
):
    """Apply NanoVNA V2 enhanced-response T/R calibration to one sweep.

    The sequence mirrors the on-screen firmware algorithm:

    1. subtract DUT leakage using the DUT's raw S11;
    2. SOL-calibrate DUT S11;
    3. subtract THRU leakage using the THRU's *own* raw S11;
    4. SOL-calibrate the stored THRU S11;
    5. calculate source-match loop gains for DUT and THRU;
    6. scale the THRU reference by ``dut_gain / thru_gain``;
    7. divide corrected DUT transmission by that enhanced THRU reference.

    Algebraically the final step is::

        dut_t  = raw21  - leakage(raw11)
        thru_t = thru21 - leakage(thru11)

        G_dut  = 1 / (1 - e11*S11_dut)
        G_thru = 1 / (1 - e11*S11_thru)

        S21 = dut_t / (thru_t * G_dut/G_thru)

    This corrects the Port-1 source-match interaction that becomes large when
    |S11| is close to one.  It is still not a full bidirectional two-port
    calibration, so a Port-2 pad remains useful for load-match suppression.
    """
    # Calibrated DUT reflection.
    s11 = _apply_reflection_cal(raw11, e00, e11, etrack)

    # Correct DUT and THRU transmission for leakage using their respective raw
    # reflection states.
    dut_leak = leak0 + raw11 * leak_r
    thru_leak = leak0 + thru11 * leak_r
    dut_t = raw21 - dut_leak
    thru_t = thru21 - thru_leak

    if np.any(np.abs(thru_t) < 1e-15):
        raise RuntimeError("Singular THRU transmission reference")

    # The THRU itself is not assumed to have exactly zero S11.  Its measured
    # reflection is corrected with the same SOL terms and used as the reference
    # state in the enhanced-response source-match correction.
    thru_s11 = _apply_reflection_cal(thru11, e00, e11, etrack)
    dut_gain = _source_match_gain(s11, e11)
    thru_gain = _source_match_gain(thru_s11, e11)

    enhanced_ref = thru_t * dut_gain / thru_gain
    if np.any(np.abs(enhanced_ref) < 1e-15):
        raise RuntimeError("Singular enhanced THRU reference")

    s21 = dut_t / enhanced_ref
    return s11, s21


def save_cal(path, f, arr, start, stop, points, info):
    """Save calibration coefficient arrays plus a JSON metadata sidecar."""
    p = Path(path)
    if p.suffix != ".npz":
        p = Path(str(p) + ".npz")

    # NPZ preserves complex numpy arrays without precision loss or a custom
    # serialization format.
    np.savez_compressed(
        p,
        frequency_hz=f,
        e00=arr[0],
        e11=arr[1],
        etrack=arr[2],
        leak0=arr[3],
        leak_r=arr[4],
        thru11=arr[5],
        thru21=arr[6],
    )

    # The JSON file is deliberately redundant: it is intended for a person to
    # inspect and records which sweep/grid/device the coefficient file belongs
    # to.
    meta = {
        "format": "unun-vna-saa2-enhanced-response",
        "version": 1,
        "algorithm": "nanovna-v2-enhanced-response",
        "start_hz": start,
        "stop_hz": stop,
        "points": points,
        "standards": ["short", "open", "load", "thru"],
        "device_info": {str(k): str(v) for k, v in info.items()},
    }
    Path(str(p) + ".json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    return p


def load_cal(path):
    """Load the enhanced-response calibration arrays used by ``apply_cal``."""
    z = np.load(path)
    required = (
        "frequency_hz",
        "e00",
        "e11",
        "etrack",
        "leak0",
        "leak_r",
        "thru11",
        "thru21",
    )
    missing = [k for k in required if k not in z.files]
    if missing:
        raise RuntimeError(
            "Calibration file is incomplete; missing: " + ", ".join(missing)
        )
    return {k: z[k] for k in required}


def same_grid(a, b):
    """Return True when two sweep frequency arrays are effectively identical.

    The 0.5 Hz absolute tolerance only accommodates representation/rounding
    details.  This version intentionally avoids interpolation of calibration
    coefficients; measurement and calibration therefore use the same sweep.
    """
    return len(a) == len(b) and np.allclose(a, b, rtol=0, atol=0.5)


def acquire_calibrated(v, cal_path):
    """Acquire one sweep on exactly the grid stored in ``cal_path``."""
    cal = load_cal(cal_path)
    cf = cal["frequency_hz"]

    # Force the measurement onto the calibration grid.  This avoids any hidden
    # interpolation of complex error terms.
    f, r11, r21 = sweep_raw(v, float(cf[0]), float(cf[-1]), len(cf))
    if not same_grid(f, cf):
        raise RuntimeError("Sweep grid differs from calibration grid")

    s11, s21 = apply_cal(
        r11,
        r21,
        cal["e00"],
        cal["e11"],
        cal["etrack"],
        cal["leak0"],
        cal["leak_r"],
        cal["thru11"],
        cal["thru21"],
    )
    return f, s11, s21


def save_measurement(path, f, s11, s21, kind, extra=None):
    """Save calibrated measurement data and a small metadata sidecar."""
    p = Path(path)
    if p.suffix != ".npz":
        p = Path(str(p) + ".npz")

    np.savez_compressed(p, frequency_hz=f, s11=s11, s21=s21)

    meta = {"format": "unun-vna-measurement", "version": 1, "kind": kind}
    if extra:
        meta.update(extra)
    Path(str(p) + ".json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    return p


def load_measurement(path):
    """Load calibrated S11/S21 data saved by :func:`save_measurement`."""
    z = np.load(path)
    return z["frequency_hz"], z["s11"], z["s21"]


def interp_real(f, values, target):
    """Linear interpolation used only for human-readable report frequencies."""
    return np.interp(target, f, values)


def interp_complex(f, values, target):
    """Interpolate real/imaginary parts separately for report output.

    The stored full-resolution data are *not* modified; interpolation is only
    used to print values at convenient band-center frequencies.
    """
    return np.interp(target, f, values.real) + 1j * np.interp(
        target, f, values.imag
    )


def cmd_info(a):
    """Implementation of ``unun-vna info``."""
    v = None
    try:
        v = open_vna(a.index)
        for k, val in v.info().items():
            print(f"{k}: {val}")
        return 0
    finally:
        if v is not None:
            try:
                v.kill()
            except Exception:
                pass


def cmd_cal(a):
    """Interactive SHORT/OPEN/LOAD/THRU calibration.

    The physical reference planes are the ends of whatever cables/adapters are
    left connected during calibration.  If a 10 dB pad is used at Port 2, it
    must already be installed and must remain there for all later measurements
    made with this calibration.
    """
    v = None
    try:
        v = open_vna(a.index)
        info = v.info()

        print("Device:")
        for k, val in info.items():
            print(f"  {k}: {val}")
        print(
            f"\nS-A-A-2 T/R SOLT: "
            f"{a.start/1e6:g}..{a.stop/1e6:g} MHz, {a.points} points"
        )
        print("Enhanced-response T/R calibration; no separate ISOLATION standard is used.")
        print("SHORT and OPEN S21 are used to derive the leakage model.")

        # Only Port 1 needs the three reflection standards because the S-A-A-2
        # does not perform a reverse/full-2-port calibration.
        fs, short11, short21 = capture(
            v,
            "1/4 SHORT — Port 1",
            "Connect SHORT at the Port-1 calibration plane.",
            a.start,
            a.stop,
            a.points,
        )
        fo, open11, open21 = capture(
            v,
            "2/4 OPEN — Port 1",
            "Connect OPEN at the Port-1 calibration plane.",
            a.start,
            a.stop,
            a.points,
        )
        fl, load11, _ = capture(
            v,
            "3/4 LOAD — Port 1",
            "Connect 50-ohm LOAD at the Port-1 calibration plane.",
            a.start,
            a.stop,
            a.points,
        )

        # THRU captures both reflection and transmission of the complete
        # forward path, including the fixed Port-2 attenuator if present.  THRU
        # S11 is required by the enhanced-response source-match correction.
        ft, thru11, thru21 = capture(
            v,
            "4/4 THRU",
            "Connect Port 1 and Port 2 calibration planes directly.",
            a.start,
            a.stop,
            a.points,
        )

        # Every coefficient is frequency-specific; mixing grids would silently
        # corrupt the calibration.
        for x in (fo, fl, ft):
            if not same_grid(fs, x):
                raise RuntimeError(
                    "Calibration sweeps returned different frequency grids"
                )

        arr = compute_cal(
            short11,
            short21,
            open11,
            open21,
            load11,
            thru11,
            thru21,
        )
        p = save_cal(
            a.output, fs, arr, a.start, a.stop, a.points, info
        )
        print(f"\nCalibration saved: {p}")
        print(f"Metadata saved:    {p}.json")

        # Immediately re-measure the still-connected THRU.  This is a necessary
        # sanity check of the calibration.  A highly reflecting series fixture
        # is still the more demanding validation because it exercises the
        # enhanced-response source-match correction and residual Port-2 match.
        f, r11, r21 = sweep_raw(v, a.start, a.stop, a.points)
        c11, c21 = apply_cal(r11, r21, *arr)
        d21 = np.array([db20(x) for x in c21])
        d11 = np.array([db20(x) for x in c11])
        vf = Path(str(p) + ".through.csv")
        write_s_csv(vf, f, c11, c21)

        print("\nTHRU verification:")
        print(
            f"  S21 min/mean/max: {d21.min():+.4f} / "
            f"{d21.mean():+.4f} / {d21.max():+.4f} dB"
        )
        print(
            f"  S11 median/worst: {np.median(d11):+.2f} / "
            f"{d11.max():+.2f} dB"
        )
        print(f"  CSV: {vf}")
        return 0
    finally:
        if v is not None:
            try:
                v.kill()
            except Exception:
                pass


def cmd_sweep(a):
    """General-purpose calibrated or uncalibrated S11/S21 CSV sweep."""
    v = None
    try:
        v = open_vna(a.index)

        if a.cal:
            cal = load_cal(a.cal)
            cf = cal["frequency_hz"]

            # When the user omits sweep settings, inherit the exact calibration
            # grid.  Explicit settings are allowed only if they still match.
            start = float(cf[0]) if a.start is None else a.start
            stop = float(cf[-1]) if a.stop is None else a.stop
            points = len(cf) if a.points is None else a.points
        else:
            cal = None
            start = DEFAULT_START if a.start is None else a.start
            stop = DEFAULT_STOP if a.stop is None else a.stop
            points = DEFAULT_POINTS if a.points is None else a.points

        f, r11, r21 = sweep_raw(v, start, stop, points)

        if cal:
            if not same_grid(f, cal["frequency_hz"]):
                raise RuntimeError("Use the same sweep grid as the calibration")
            s11, s21 = apply_cal(
                r11,
                r21,
                cal["e00"],
                cal["e11"],
                cal["etrack"],
                cal["leak0"],
                cal["leak_r"],
                cal["thru11"],
                cal["thru21"],
            )
            print(f"Applied S-A-A-2 T/R SOLT calibration: {a.cal}")
        else:
            s11, s21 = r11, r21
            print("UNCALIBRATED sweep")

        write_s_csv(a.output, f, s11, s21)

        # Print a small preview rather than flooding the terminal with all 501
        # points.  The complete data are always in the CSV.
        idx = sorted(
            set([0, 1, 2, len(f) // 2, len(f) - 3, len(f) - 2, len(f) - 1])
        )
        for i in idx:
            print(
                f"{f[i]/1e6:9.4f} MHz  "
                f"S11={s11[i].real:+.6e}{s11[i].imag:+.6e}j "
                f"({db20(s11[i]):+.3f} dB)  "
                f"S21={s21[i].real:+.6e}{s21[i].imag:+.6e}j "
                f"({db20(s21[i]):+.3f} dB)"
            )
        print(f"CSV written to: {a.output}")
        return 0
    finally:
        if v is not None:
            try:
                v.kill()
            except Exception:
                pass


def cmd_measure(a):
    """Acquire the complete UnUn + resistor-fixture DUT measurement.

    This command does not calculate efficiency immediately.  It preserves the
    calibrated complex S-parameters so the analysis can later be rerun with a
    revised fixture model without repeating the RF measurement.
    """
    v = None
    try:
        v = open_vna(a.index)
        f, s11, s21 = acquire_calibrated(v, a.cal)

        p = save_measurement(
            a.output,
            f,
            s11,
            s21,
            "dut",
            {"calibration": str(a.cal)},
        )

        csv_path = Path(str(p) + ".csv")
        write_s_csv(csv_path, f, s11, s21)

        print(f"Measurement saved: {p}")
        print(f"CSV:               {csv_path}")
        return 0
    finally:
        if v is not None:
            try:
                v.kill()
            except Exception:
                pass


def cmd_fixture(a):
    """Measure and validate an arbitrary series-load fixture.

    This step is deliberately separate from DUT measurement.  No nominal
    resistance is assumed; the fixture's complex series impedance is extracted
    from calibrated S11/S21 at every frequency.

    The two independent extraction formulae and the identity S11+S21=1 are used
    as diagnostics.  If those checks fail badly, the fixture cannot safely be
    represented by the simple series-impedance model used by the efficiency
    calculation.
    """
    v = None
    try:
        v = open_vna(a.index)
        f, s11, s21 = acquire_calibrated(v, a.cal)

        # Mathematical extraction/validation lives in analysis.py.
        model = series_fixture_impedance(s11, s21, a.z0)
        z = model["z_series"]

        p = Path(a.output)
        if p.suffix != ".npz":
            p = Path(str(p) + ".npz")

        # Store both raw calibrated S-parameters and all derived diagnostics.
        # This makes the fixture file auditable and allows the extraction method
        # to be revisited later without re-measuring the hardware.
        np.savez_compressed(
            p,
            frequency_hz=f,
            s11=s11,
            s21=s21,
            z_series=z,
            z_from_ratio=model["z_from_ratio"],
            z_from_s21=model["z_from_s21"],
            identity_error=model["identity_error"],
            model_error=model["model_error"],
            z0=np.array(a.z0),
        )

        Path(str(p) + ".json").write_text(
            json.dumps(
                {
                    "format": "unun-vna-series-fixture",
                    "version": 1,
                    "calibration": str(a.cal),
                    "z0_ohm": a.z0,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        # The CSV is intended for plotting/inspection in external tools.
        csv_path = Path(str(p) + ".csv")
        with csv_path.open("w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(
                [
                    "frequency_hz",
                    "s11_db",
                    "s21_db",
                    "z_series_real_ohm",
                    "z_series_imag_ohm",
                    "z_ratio_real_ohm",
                    "z_ratio_imag_ohm",
                    "z_s21_real_ohm",
                    "z_s21_imag_ohm",
                    "series_identity_error_abs",
                    "model_error_abs",
                ]
            )
            for i in range(len(f)):
                w.writerow(
                    [
                        int(round(f[i])),
                        db20(s11[i]),
                        db20(s21[i]),
                        z[i].real,
                        z[i].imag,
                        model["z_from_ratio"][i].real,
                        model["z_from_ratio"][i].imag,
                        model["z_from_s21"][i].real,
                        model["z_from_s21"][i].imag,
                        abs(model["identity_error"][i]),
                        model["model_error"][i],
                    ]
                )

        # Summarize how closely the physical PCB behaves like the intended
        # single series element.  Small residuals are more important here than
        # agreement with any nominal resistance.
        ident = np.abs(model["identity_error"])
        disagreement = np.abs(
            model["z_from_ratio"] - model["z_from_s21"]
        )

        print(f"Fixture saved: {p}")
        print(f"CSV:           {csv_path}")
        print("\nPure-series model diagnostics:")
        print(f"  Re(Zs) median: {np.median(z.real):.2f} ohm")
        print(f"  Im(Zs) median: {np.median(z.imag):+.2f} ohm")
        print(
            f"  |S11+S21-1| median/max: "
            f"{np.median(ident):.4g} / {np.max(ident):.4g}"
        )
        print(
            f"  two Z estimators disagreement median/max: "
            f"{np.median(disagreement):.2f} / "
            f"{np.max(disagreement):.2f} ohm"
        )

        # Ratio affects diagnostics only. The measured fixture impedance remains
        # authoritative regardless of the configured transformer ratio.
        ref = ideal_transformer_reference(z, a.ratio, a.z0)
        print(
            f"\nTransformer reference: low:high impedance ratio "
            f"1:{a.ratio:g}  (turns ratio 1:{ref['turns_ratio']:.5g})"
        )
        print(
            f"  target total high-side resistance: "
            f"{ref['target_high_side_resistance']:.2f} ohm"
        )
        print(
            f"  target pure-series fixture resistance: "
            f"{ref['target_series_resistance']:.2f} ohm"
        )
        print("  NOTE: target values are diagnostics only; measured Zs is authoritative.")

        for target in parse_freq_list(a.bands):
            v11 = interp_complex(f, s11, target)
            v21 = interp_complex(f, s21, target)
            zz = interp_complex(f, z, target)
            ident_t = interp_real(f, ident, target)
            print(
                f"  {target/1e6:6.3f} MHz: "
                f"S11={db20(v11):+.3f} dB, "
                f"S21={db20(v21):+.3f} dB, "
                f"Zs={zz.real:.2f}{zz.imag:+.2f}j ohm, "
                f"series-error={ident_t:.4g}"
            )

        if np.max(ident) > a.warn_series_error:
            print(
                "\nWARNING: Fixture departs noticeably from a pure series "
                "impedance at some frequencies. Efficiency results should be "
                "treated cautiously there.",
                file=sys.stderr,
            )
        return 0
    finally:
        if v is not None:
            try:
                v.kill()
            except Exception:
                pass


def cmd_analyze(a):
    """Combine DUT and fixture measurements and calculate UnUn efficiency."""
    # Load the complete calibrated DUT measurement.
    f_dut, s11, s21 = load_measurement(a.dut)

    # Load the independently characterized series-load fixture.
    fixture = np.load(a.fixture)
    f_fix = fixture["frequency_hz"]

    # Point-for-point arithmetic is intentional.  We reject mismatched grids
    # rather than silently interpolating measurement data.
    if not same_grid(f_dut, f_fix):
        raise RuntimeError(
            "DUT and fixture must use the same calibrated frequency grid"
        )

    # The fixture command stores the extracted complex series impedance
    # explicitly; analysis consumes exactly that measured model.
    z_series = fixture["z_series"]

    # The actual power calculation is kept in analysis.py so it can be tested
    # independently of CLI/file I/O.
    result = unun_efficiency(s11, s21, z_series, a.z0)
    eta = result["efficiency"]
    loss_db = result["loss_db"]

    # Ratio-dependent values are reference diagnostics only. The measured
    # efficiency above intentionally does not use a.ratio.
    ideal_ref = ideal_transformer_reference(z_series, a.ratio, a.z0)
    ideal_s11 = ideal_ref["gamma_in"]
    ideal_s21_mag = ideal_ref["s21_magnitude"]

    # VSWR is a direct re-expression of the calibrated Port-1 reflection.
    # The ideal-reference VSWR shows the load mismatch that a lossless
    # transformer of the configured ratio would have with the actually
    # measured fixture.
    vswr = vswr_from_s11(s11)
    ideal_vswr = vswr_from_s11(ideal_s11)

    # Write the full frequency-resolved result.  The complex fixture impedance
    # is included so every efficiency point can be audited later.
    out = Path(a.output)
    with out.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "frequency_hz",
                "s11_db",
                "vswr",
                "s21_db",
                "fixture_r_ohm",
                "fixture_x_ohm",
                "impedance_ratio_high_over_low",
                "target_high_side_resistance_ohm",
                "target_series_resistance_ohm",
                "ideal_lossless_s11_db_for_measured_fixture",
                "ideal_lossless_vswr_for_measured_fixture",
                "ideal_lossless_s21_db_for_measured_fixture",
                "accepted_fraction",
                "delivered_fraction",
                "efficiency",
                "efficiency_percent",
                "transformer_loss_db",
            ]
        )
        for i in range(len(f_dut)):
            w.writerow(
                [
                    int(round(f_dut[i])),
                    db20(s11[i]),
                    vswr[i],
                    db20(s21[i]),
                    z_series[i].real,
                    z_series[i].imag,
                    a.ratio,
                    ideal_ref["target_high_side_resistance"],
                    ideal_ref["target_series_resistance"],
                    db20(ideal_s11[i]),
                    ideal_vswr[i],
                    db20(ideal_s21_mag[i]),
                    result["accepted_fraction"][i],
                    result["delivered_fraction"][i],
                    eta[i],
                    eta[i] * 100.0,
                    loss_db[i],
                ]
            )

    print(f"Analysis CSV: {out}")
    print(
        f"Transformer ratio: low:high impedance = 1:{a.ratio:g} "
        f"(turns ratio 1:{ideal_ref['turns_ratio']:.5g})"
    )
    print(
        f"Ideal purely-resistive matching load: "
        f"total {ideal_ref['target_high_side_resistance']:.2f} ohm, "
        f"series fixture {ideal_ref['target_series_resistance']:.2f} ohm"
    )
    print("Measured fixture impedance is used for all efficiency calculations.")

    print("\nSelected frequencies:")
    print(
        " frequency     S11      VSWR      S21       ideal S11/VSWR/S21"
        "            fixture Zs             efficiency   loss"
    )

    # Interpolate only for compact reporting at convenient band frequencies.
    # The CSV itself contains the original grid without interpolation.
    for target in parse_freq_list(a.bands):
        v11 = interp_complex(f_dut, s11, target)
        v21 = interp_complex(f_dut, s21, target)
        zz = interp_complex(f_dut, z_series, target)
        ee = interp_real(f_dut, eta, target)
        ll = interp_real(f_dut, loss_db, target)
        vv = interp_real(f_dut, vswr, target)
        i11 = interp_complex(f_dut, ideal_s11, target)
        iv = interp_real(f_dut, ideal_vswr, target)
        i21 = interp_real(f_dut, ideal_s21_mag, target)
        print(
            f" {target/1e6:7.3f} MHz  "
            f"{db20(v11):7.2f} dB  "
            f"{vv:7.3f}  "
            f"{db20(v21):7.2f} dB  "
            f"{db20(i11):7.2f}/{iv:6.3f}/{db20(i21):7.2f} dB  "
            f"{zz.real:7.1f}{zz.imag:+7.1f}j ohm  "
            f"{ee*100:7.2f} %   {ll:6.3f} dB"
        )

    # eta > 100% is a useful metrology diagnostic.  A tiny amount may result
    # from residual noise/calibration error, but a substantial excess indicates
    # that load-match error or the fixture model is still inadequate.
    finite = np.isfinite(eta)
    if np.any(finite & (eta > 1.02)):
        print(
            "\nWARNING: Efficiency exceeds 102% at some frequencies. "
            "This indicates residual calibration/load-match error or a fixture "
            "that is not adequately represented by the series model.",
            file=sys.stderr,
        )
    return 0


def make_parser():
    """Construct the command-line interface."""
    p = argparse.ArgumentParser(
        prog="unun-vna",
        description=(
            "S-A-A-2/NanoVNA V2 tool for one-transformer UnUn efficiency "
            "measurements"
        ),
    )
    s = p.add_subparsers(dest="cmd", required=True)

    q = s.add_parser("info", help="Show detected VNA information")
    q.add_argument("--index", type=int, default=0)

    q = s.add_parser(
        "calibrate",
        help="Perform host-side SHORT/OPEN/LOAD/THRU T/R calibration",
    )
    q.add_argument("--index", type=int, default=0)
    q.add_argument("--start", type=parse_freq, default=DEFAULT_START)
    q.add_argument("--stop", type=parse_freq, default=DEFAULT_STOP)
    q.add_argument("--points", type=int, default=DEFAULT_POINTS)
    q.add_argument("-o", "--output", default="hf35")

    q = s.add_parser(
        "sweep",
        help="Acquire a generic calibrated or uncalibrated S11/S21 sweep",
    )
    q.add_argument("--index", type=int, default=0)
    q.add_argument("--start", type=parse_freq)
    q.add_argument("--stop", type=parse_freq)
    q.add_argument("--points", type=int)
    q.add_argument("--cal", type=Path)
    q.add_argument("-o", "--output", default="saa2_sweep.csv")

    q = s.add_parser(
        "fixture",
        help="Measure and characterize an arbitrary series-load fixture",
    )
    q.add_argument("--index", type=int, default=0)
    q.add_argument("--cal", type=Path, required=True)
    q.add_argument("-o", "--output", default="fixture")
    q.add_argument("--z0", type=float, default=DEFAULT_Z0)
    q.add_argument(
        "--ratio",
        type=parse_impedance_ratio,
        default=DEFAULT_RATIO,
        help=(
            "transformer impedance ratio high/low, e.g. 49 or 1:49 "
            "(default: 49); diagnostics only"
        ),
    )
    q.add_argument("--bands", default=DEFAULT_BANDS)
    q.add_argument(
        "--warn-series-error",
        type=float,
        default=0.01,
        help="Warn if max |S11+S21-1| exceeds this value",
    )

    q = s.add_parser(
        "measure",
        help="Acquire calibrated UnUn + fixture data for later analysis",
    )
    q.add_argument("--index", type=int, default=0)
    q.add_argument("--cal", type=Path, required=True)
    q.add_argument("-o", "--output", default="dut")

    q = s.add_parser(
        "analyze",
        help="Calculate UnUn efficiency using measured fixture impedance",
    )
    q.add_argument("--dut", type=Path, required=True)
    q.add_argument("--fixture", type=Path, required=True)
    q.add_argument("--z0", type=float, default=DEFAULT_Z0)
    q.add_argument(
        "--ratio",
        type=parse_impedance_ratio,
        default=DEFAULT_RATIO,
        help=(
            "transformer impedance ratio high/low, e.g. 49 or 1:49 "
            "(default: 49)"
        ),
    )
    q.add_argument("--bands", default=DEFAULT_BANDS)
    q.add_argument("-o", "--output", default="unun-efficiency.csv")

    return p


def main():
    """CLI entry point installed as the ``unun-vna`` console command."""
    a = make_parser().parse_args()
    try:
        if a.cmd == "info":
            return cmd_info(a)
        if a.cmd == "calibrate":
            return cmd_cal(a)
        if a.cmd == "sweep":
            return cmd_sweep(a)
        if a.cmd == "fixture":
            return cmd_fixture(a)
        if a.cmd == "measure":
            return cmd_measure(a)
        if a.cmd == "analyze":
            return cmd_analyze(a)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
