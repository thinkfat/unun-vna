"""Visualization of ``unun-vna analyze`` CSV output.

The plots are deliberately aimed at the practical questions an amateur-radio
operator usually asks first:

* How well is the transformer matched?  -> VSWR
* How much accepted RF power reaches the load? -> transformer efficiency

The plotting layer consumes only the CSV written by ``analyze``.  It performs
no RF calculations and therefore cannot change the measurement result.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


# Overall amateur allocations used for the visual band markers.  These are the
# HF allocations relevant to the default 1..35 MHz sweep in IARU Region 1.
IARU_REGION1_BANDS = (
    ("40 m", 7.000, 7.200),
    ("20 m", 14.000, 14.350),
    ("15 m", 21.000, 21.450),
    ("10 m", 28.000, 29.700),
)

SUPPORTED_OUTPUT_SUFFIXES = {".png", ".svg", ".pdf"}


@dataclass(frozen=True)
class AnalysisSeries:
    """One analyzed transformer sweep ready for plotting."""

    path: Path
    label: str
    frequency_mhz: np.ndarray
    vswr: np.ndarray
    efficiency_percent: np.ndarray


def _default_label(path: Path) -> str:
    """Create a compact legend label from an analysis CSV filename."""
    stem = path.stem
    for suffix in ("-efficiency", "_efficiency", "-analysis", "_analysis"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem or path.name


def read_analysis_csv(path: Path | str, label: str | None = None) -> AnalysisSeries:
    """Read the columns required for the user-facing performance plot.

    Required columns are ``frequency_hz`` and ``vswr`` plus either
    ``efficiency_percent`` or the unit-fraction column ``efficiency``.
    Additional columns written by ``analyze`` are intentionally ignored.
    """
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"analysis CSV does not exist: {p}")

    frequency: list[float] = []
    vswr: list[float] = []
    efficiency: list[float] = []

    with p.open("r", newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        fields = set(reader.fieldnames or [])
        required = {"frequency_hz", "vswr"}
        missing = required - fields
        if missing:
            raise ValueError(
                f"{p}: not an analyze CSV; missing column(s): "
                + ", ".join(sorted(missing))
            )

        if "efficiency_percent" in fields:
            efficiency_key = "efficiency_percent"
            efficiency_scale = 1.0
        elif "efficiency" in fields:
            efficiency_key = "efficiency"
            efficiency_scale = 100.0
        else:
            raise ValueError(
                f"{p}: missing efficiency_percent (or efficiency) column"
            )

        for line_no, row in enumerate(reader, start=2):
            try:
                f_hz = float(row["frequency_hz"])
                v = float(row["vswr"])
                e = float(row[efficiency_key]) * efficiency_scale
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{p}:{line_no}: invalid numeric data") from exc

            if not math.isfinite(f_hz):
                continue
            frequency.append(f_hz / 1e6)
            vswr.append(v)
            efficiency.append(e)

    if len(frequency) < 2:
        raise ValueError(f"{p}: fewer than two valid frequency points")

    f = np.asarray(frequency, dtype=float)
    v = np.asarray(vswr, dtype=float)
    e = np.asarray(efficiency, dtype=float)

    # Analysis normally writes an ascending sweep.  Sorting here makes the
    # standalone plot command robust to externally reordered CSV rows.
    order = np.argsort(f)
    return AnalysisSeries(
        path=p,
        label=label if label is not None else _default_label(p),
        frequency_mhz=f[order],
        vswr=v[order],
        efficiency_percent=e[order],
    )


def _shade_region1_bands(ax, x_min: float, x_max: float, add_labels: bool) -> None:
    """Mark the 40/20/15/10 m IARU Region-1 allocations on one axis."""
    for name, start, stop in IARU_REGION1_BANDS:
        if stop < x_min or start > x_max:
            continue
        ax.axvspan(max(start, x_min), min(stop, x_max), color="0.92", zorder=0)
        if add_labels:
            mid = (max(start, x_min) + min(stop, x_max)) / 2.0
            # The lower-band allocations are very narrow on a 1..35 MHz plot;
            # vertical labels stay readable without artificially widening them.
            ax.text(
                mid,
                0.04,
                name,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=8,
                color="0.35",
                transform=ax.get_xaxis_transform(),
            )


def _automatic_efficiency_min(series: Sequence[AnalysisSeries]) -> float:
    """Choose a non-misleading rounded lower limit while retaining all data."""
    finite_chunks = [
        s.efficiency_percent[np.isfinite(s.efficiency_percent)] for s in series
    ]
    finite_chunks = [x for x in finite_chunks if len(x)]
    if not finite_chunks:
        return 0.0
    minimum = float(np.min(np.concatenate(finite_chunks)))
    # Round down to a 10-percentage-point boundary.  Never use a truncated
    # baseline above zero if the data themselves reach zero or below.
    return max(0.0, min(90.0, math.floor(minimum / 10.0) * 10.0))


def plot_analysis_files(
    inputs: Sequence[Path | str],
    output: Path | str,
    *,
    labels: Sequence[str] | None = None,
    title: str = "UnUn performance",
    show_bands: bool = True,
    vswr_max: float = 3.0,
    efficiency_min: float | None = None,
    dpi: int = 160,
) -> Path:
    """Create the standard two-panel VSWR/efficiency visualization.

    ``inputs`` may contain multiple analysis CSVs.  Series are plotted in the
    same order in both panels so Matplotlib's color cycle assigns matching
    colors without the code hard-coding a palette.
    """
    if not inputs:
        raise ValueError("at least one analysis CSV is required")
    if vswr_max <= 1.0 or not math.isfinite(vswr_max):
        raise ValueError("vswr_max must be finite and greater than 1")
    if efficiency_min is not None and not (0.0 <= efficiency_min < 100.0):
        raise ValueError("efficiency_min must be in the range 0 <= value < 100")
    if dpi <= 0:
        raise ValueError("dpi must be > 0")

    if labels is not None and len(labels) not in (0, len(inputs)):
        raise ValueError(
            "when --label is used, provide exactly one label for each input CSV"
        )

    resolved_labels = list(labels or [])
    series = [
        read_analysis_csv(path, resolved_labels[i] if resolved_labels else None)
        for i, path in enumerate(inputs)
    ]

    out = Path(output)
    suffix = out.suffix.lower()
    if suffix not in SUPPORTED_OUTPUT_SUFFIXES:
        raise ValueError(
            "plot output must end in .png, .svg or .pdf "
            f"(got {out.name!r})"
        )
    out.parent.mkdir(parents=True, exist_ok=True)

    # Force a non-interactive backend so the CLI works over SSH and on systems
    # without a desktop session.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_vswr, ax_eff) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(10.5, 7.2),
        constrained_layout=True,
    )

    x_min = min(float(np.min(s.frequency_mhz)) for s in series)
    x_max = max(float(np.max(s.frequency_mhz)) for s in series)

    if show_bands:
        _shade_region1_bands(ax_vswr, x_min, x_max, add_labels=True)
        _shade_region1_bands(ax_eff, x_min, x_max, add_labels=False)

    for s in series:
        ax_vswr.plot(s.frequency_mhz, s.vswr, label=s.label, linewidth=1.6)
        ax_eff.plot(
            s.frequency_mhz,
            s.efficiency_percent,
            label=s.label,
            linewidth=1.6,
        )

    # VSWR is intentionally capped visually: a single badly mismatched point
    # must not make the useful 1:1..2:1 region unreadable.  The underlying CSV
    # retains the full value and users can override the limit with --vswr-max.
    ax_vswr.set_ylim(1.0, vswr_max)
    ax_vswr.set_ylabel("VSWR")
    ax_vswr.axhline(1.5, linestyle="--", linewidth=0.9, color="0.45")
    ax_vswr.axhline(2.0, linestyle=":", linewidth=0.9, color="0.45")
    ax_vswr.grid(True, alpha=0.25)

    eff_min = (
        _automatic_efficiency_min(series)
        if efficiency_min is None
        else float(efficiency_min)
    )
    ax_eff.set_ylim(eff_min, 100.0)
    ax_eff.set_ylabel("Efficiency (%)")
    ax_eff.set_xlabel("Frequency (MHz)")
    ax_eff.grid(True, alpha=0.25)

    ax_eff.set_xlim(x_min, x_max)

    # A legend is useful even for one series because plots are frequently
    # copied out of their original measurement-directory context.
    ax_vswr.legend(loc="best")
    ax_eff.legend(loc="best")

    fig.suptitle(title)
    fig.savefig(out, dpi=dpi if suffix == ".png" else None)
    plt.close(fig)
    return out
