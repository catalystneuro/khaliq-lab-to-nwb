"""
Utilities to detect saturated stretches in intracellular recordings stored in NWB files and visualize them.

Usage
-----
Run with uv (recommended in this repository):

    UV_CACHE_DIR=.uv_cache uv run python -m khaliq_lab_to_nwb.conversion_abf.detect_saturation \\
        nwbfiles/subject_1.nwb \\
        --output-dir build/plots

The script scans every acquisition `TimeSeries`, flags saturated intervals, and writes an SVG plot with all
traces overlaid and translucent boxes around the regions that were flagged as corrupted.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBHDF5IO
from pynwb.base import TimeSeries


@dataclass
class DetectionResult:
    series_name: str
    times: np.ndarray
    values: np.ndarray
    unit: str
    rectangles: list[tuple[float, float, float, float]]


def _take_every_n(arr: np.ndarray, n: int) -> np.ndarray:
    """Return evenly spaced indices to downsample an array without losing endpoints."""
    if len(arr) <= n:
        return np.arange(len(arr))
    step = (len(arr) - 1) / (n - 1)
    return np.round(np.arange(n) * step).astype(int)


def _load_series(
    series: TimeSeries,
    max_points: int = 8000,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return relative timestamps, values, and absolute start time for plotting."""
    data = np.array(series.data[:], dtype=np.float64)
    if series.timestamps is not None:
        times = np.array(series.timestamps[:], dtype=np.float64)
    else:
        rate = float(series.rate)
        start = float(series.starting_time)
        times = start + np.arange(data.size, dtype=np.float64) / rate

    if times.size == 0:
        return times, data, 0.0

    start_time = float(times[0])
    rel_times = times - start_time

    if data.size > max_points:
        idx = _take_every_n(data, max_points)
        return rel_times[idx], data[idx], start_time
    return rel_times, data, start_time


def detect_saturation_segments(
    data: np.ndarray,
    rate: float,
    start_time: float,
    amp_factor: float = 0.2,
) -> list[tuple[float, float, float, float]]:
    """
    Detect saturated stretches based on absolute voltage thresholds.

    Normal intracellular recordings should be within ±200 mV. Values beyond this
    indicate amplifier saturation or data corruption.

    Segments are identified by grouping voltage exceedances that occur within 100ms
    of each other. Single-point spikes are filtered out.

    Parameters
    ----------
    data : np.ndarray
        Voltage samples (volts).
    rate : float
        Sampling rate in Hz.
    start_time : float
        Starting time offset of the series in seconds.
    amp_factor : float
        Absolute voltage threshold in volts (default 0.2 V = 200 mV).
        Values exceeding ±amp_factor are considered saturated.

    Returns
    -------
    list[tuple[float, float, float, float]]
        List of saturation segments, each as (start_time, stop_time, voltage_min, voltage_max).
    """
    if data.size == 0:
        return []

    # Use absolute voltage threshold: ±200 mV (0.2 V)
    # amp_factor is repurposed as the absolute threshold in volts
    threshold_V = amp_factor if amp_factor < 1.0 else 0.2

    # Find points that exceed the threshold
    extreme_mask = (data > threshold_V) | (data < -threshold_V)

    if not np.any(extreme_mask):
        return []

    # Find continuous segments where voltage exceeds threshold
    extreme_indices = np.where(extreme_mask)[0]

    if len(extreme_indices) == 0:
        return []

    # Group exceedances into segments based on temporal proximity
    # Points within 100ms of each other belong to the same corrupted segment
    segments: list[list[int]] = []
    current_segment = [extreme_indices[0]]
    gap_threshold_samples = int(0.1 * rate)  # 100ms in samples

    for idx in extreme_indices[1:]:
        if idx - current_segment[-1] <= gap_threshold_samples:
            # This point is within 100ms of the previous one - same segment
            current_segment.append(idx)
        else:
            # Gap > 100ms - start a new segment
            segments.append(current_segment)
            current_segment = [idx]

    # Don't forget the last segment
    segments.append(current_segment)

    # Filter out segments with only 1 point (isolated spikes)
    valid_segments = [seg for seg in segments if len(seg) > 1]

    if not valid_segments:
        # Only isolated single-point spikes - not corruption
        return []

    # Convert each valid segment to a time range
    results: list[tuple[float, float, float, float]] = []
    for segment_indices in valid_segments:
        first_idx = segment_indices[0]
        last_idx = segment_indices[-1]

        # Get voltage range for this corrupted segment
        segment_data = data[first_idx : last_idx + 1]
        seg_min = float(np.min(segment_data))
        seg_max = float(np.max(segment_data))

        # Add small epsilon if min==max
        if math.isclose(seg_min, seg_max):
            epsilon = max(0.0005 * abs(seg_min if seg_min else 1.0), 1e-6)
            seg_min -= epsilon
            seg_max += epsilon

        # Convert indices to absolute times
        time_start = start_time + first_idx / rate
        time_end = start_time + last_idx / rate

        results.append((time_start, time_end, seg_min, seg_max))

    return results


def analyze_file(
    nwb_path: Path,
    amp_factor: float,
    max_points: int,
) -> list[DetectionResult]:
    """Run saturation detection on every acquisition series in an NWB file."""
    results: list[DetectionResult] = []
    with NWBHDF5IO(nwb_path, "r") as io:
        nwbfile = io.read()
        acquisitions = list(nwbfile.acquisition.values())
        for series in acquisitions:
            if not isinstance(series, TimeSeries):
                continue

            full_data = np.array(series.data[:], dtype=np.float64)
            if series.timestamps is not None:
                start_time = float(series.timestamps[0])
            else:
                start_time = float(series.starting_time)
            rate = float(series.rate)
            rectangles = detect_saturation_segments(
                data=full_data,
                rate=rate,
                start_time=start_time,
                amp_factor=amp_factor,
            )
            rel_times, display_values, series_start = _load_series(
                series, max_points=max_points
            )
            rectangles_rel = [
                (start - series_start, stop - series_start, vmin, vmax)
                for start, stop, vmin, vmax in rectangles
            ]
            results.append(
                DetectionResult(
                    series_name=series.name,
                    times=rel_times,
                    values=display_values,
                    unit=getattr(series, "unit", "unknown"),
                    rectangles=rectangles_rel,
                )
            )
    return results


def _choose_color(index: int) -> str:
    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    return palette[index % len(palette)]


def render_png(
    results: Sequence[DetectionResult], output_path: Path, title: str
) -> None:
    """Render a PNG plot with all traces and saturation intervals marked with shaded regions using matplotlib."""
    if not results:
        raise ValueError("No acquisition TimeSeries to plot.")

    # Create figure with white background
    fig, ax = plt.subplots(figsize=(14, 8), facecolor="white")
    ax.set_facecolor("white")

    # Find the time range of corrupted segments to set plot limits
    # Show 1 second before first corruption and 1 second after last corruption
    all_corruption_times = []
    for res in results:
        for start, stop, vmin, vmax in res.rectangles:
            all_corruption_times.extend([start, stop])

    if all_corruption_times:
        first_corruption = min(all_corruption_times)
        last_corruption = max(all_corruption_times)
        time_plot_min = first_corruption - 1.0  # 1 second before
        time_plot_max = last_corruption + 1.0  # 1 second after
    else:
        # If no corruption found, use full time range
        time_plot_min = min(res.times[0] for res in results if res.times.size > 0)
        time_plot_max = max(res.times[-1] for res in results if res.times.size > 0)

    # Collect values for voltage range calculation
    all_values_mV = []
    for res in results:
        if res.times.size > 0:
            # Only include values within the plot time range
            mask = (res.times >= time_plot_min) & (res.times <= time_plot_max)
            if np.any(mask):
                all_values_mV.append(res.values[mask] * 1000)

    if all_values_mV:
        all_values_mV = np.concatenate(all_values_mV)
        value_min_mV = float(np.min(all_values_mV))
        value_max_mV = float(np.max(all_values_mV))
    else:
        value_min_mV, value_max_mV = -100, 50

    # Plot all traces in millivolts
    for idx, res in enumerate(results):
        color = _choose_color(idx)
        if res.times.size > 0:
            values_mV = res.values * 1000
            # Only plot data within the time range
            mask = (res.times >= time_plot_min) & (res.times <= time_plot_max)
            ax.plot(
                res.times[mask], values_mV[mask], color=color, linewidth=1.0, alpha=0.8
            )

    # Mark saturation intervals with shaded vertical spans
    for idx, res in enumerate(results):
        for start, stop, vmin, vmax in res.rectangles:
            ax.axvspan(
                start,
                stop,
                color="red",
                alpha=0.2,
                label="Corrupted region" if idx == 0 else "",
            )

    # Set time limits to show 1 second before and after corruption
    ax.set_xlim(time_plot_min, time_plot_max)

    # Styling with dark colors for visibility on white background
    ax.set_xlabel("Time (s)", color="black", fontsize=14)
    ax.set_ylabel("Amplitude (mV)", color="black", fontsize=14)
    ax.set_title(title, color="black", fontsize=18, pad=20)
    ax.tick_params(colors="black", labelsize=12)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#cccccc", alpha=0.5, linestyle="-", linewidth=0.5)

    # Add text annotation in mV
    ax.text(
        0.5,
        -0.12,
        f"Amplitude range: {value_min_mV:.2f} to {value_max_mV:.2f} mV",
        transform=ax.transAxes,
        ha="center",
        color="black",
        fontsize=12,
    )

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, facecolor="white", edgecolor="none")
    plt.close()


def render_svg(
    results: Sequence[DetectionResult], output_path: Path, title: str
) -> None:
    """Render an SVG overlay with all traces and saturation rectangles."""
    if not results:
        raise ValueError("No acquisition TimeSeries to plot.")

    width, height = 1400, 800
    margin_left, margin_right = 80, 40
    margin_top, margin_bottom = 60, 80
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom

    all_time_arrays = [res.times for res in results if res.times.size > 0]
    if not all_time_arrays:
        raise ValueError("No time samples found in the acquisition series.")
    all_times = np.concatenate(all_time_arrays)
    rect_time_values: list[float] = []
    for res in results:
        for start, stop, _, _ in res.rectangles:
            rect_time_values.extend([start, stop])
    if rect_time_values:
        all_times = np.concatenate(
            [all_times, np.array(rect_time_values, dtype=np.float64)]
        )

    time_min = float(np.min(all_times))
    time_max = float(np.max(all_times))
    if math.isclose(time_min, time_max):
        time_max = time_min + 1.0

    all_values = np.concatenate([res.values for res in results])
    value_min = float(np.min(all_values))
    value_max = float(np.max(all_values))

    for res in results:
        for start, stop, vmin, vmax in res.rectangles:
            value_min = min(value_min, vmin)
            value_max = max(value_max, vmax)

    if math.isclose(value_min, value_max):
        value_max = value_min + 1e-3

    def x_coord(time_value: float) -> float:
        return (
            margin_left + (time_value - time_min) / (time_max - time_min) * plot_width
        )

    def y_coord(value: float) -> float:
        return (
            margin_top
            + (1 - (value - value_min) / (value_max - value_min)) * plot_height
        )

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#111a1f"/>',
        f'<text x="{width / 2}" y="{margin_top / 2}" fill="#e5e8eb" font-size="24" text-anchor="middle">{title}</text>',
        f'<text x="{width / 2}" y="{height - margin_bottom / 3}" fill="#9fa6ad" font-size="16" text-anchor="middle">'
        f"Amplitude range: {value_min:.3f} to {value_max:.3f} ({results[0].unit})</text>",
    ]

    # Axes
    axis_color = "#4b5563"
    lines.append(
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" '
        f'y2="{margin_top + plot_height}" stroke="{axis_color}" stroke-width="1.5"/>'
    )
    lines.append(
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" '
        f'y2="{margin_top + plot_height}" stroke="{axis_color}" stroke-width="1.5"/>'
    )

    # Ticks
    tick_font = 14
    for frac in np.linspace(0.0, 1.0, 6):
        t = time_min + frac * (time_max - time_min)
        x = margin_left + frac * plot_width
        y = margin_top + plot_height
        lines.append(
            f'<line x1="{x}" y1="{y}" x2="{x}" y2="{y + 10}" stroke="{axis_color}" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{x}" y="{y + 28}" fill="#9fa6ad" font-size="{tick_font}" text-anchor="middle">{t:.2f}s</text>'
        )

    for frac in np.linspace(0.0, 1.0, 5):
        v = value_min + frac * (value_max - value_min)
        y = margin_top + (1 - frac) * plot_height
        lines.append(
            f'<line x1="{margin_left - 10}" y1="{y}" x2="{margin_left}" y2="{y}" stroke="{axis_color}" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{margin_left - 14}" y="{y + 5}" fill="#9fa6ad" font-size="{tick_font}" text-anchor="end">{v:.3f}</text>'
        )

    # Plot rectangles first (so traces render on top)
    for idx, res in enumerate(results):
        color = _choose_color(idx)
        for start, stop, vmin, vmax in res.rectangles:
            x = x_coord(start)
            y = y_coord(vmax)
            width_rect = max(x_coord(stop) - x, 2.0)
            height_rect = max(y_coord(vmin) - y, 2.0)
            lines.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{width_rect:.2f}" height="{height_rect:.2f}" '
                f'fill="{color}22" stroke="{color}" stroke-width="1.5"/>'
            )

    for idx, res in enumerate(results):
        color = _choose_color(idx)
        coords = " ".join(
            f"{x_coord(t):.2f},{y_coord(v):.2f}"
            for t, v in zip(res.times, res.values, strict=False)
        )
        lines.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="1.5" stroke-linejoin="round" '
            f'stroke-linecap="round" points="{coords}"/>'
        )

    # Legend
    legend_x = margin_left + plot_width + 10
    legend_y = margin_top + 20
    lines.append(f'<g transform="translate({legend_x}, {legend_y})">')
    lines.append(
        f'<text x="0" y="-10" fill="#d1d5db" font-size="16">Sweeps ({len(results)})</text>'
    )
    for idx, res in enumerate(
        results[:20]
    ):  # limit legend entries to first 20 for readability
        color = _choose_color(idx)
        lines.append(
            f'<rect x="0" y="{idx * 24}" width="20" height="4" fill="{color}"/>'
        )
        lines.append(
            f'<text x="26" y="{idx * 24 + 6}" fill="#cbd5f5" font-size="12">{res.series_name}</text>'
        )
    lines.append("</g>")

    lines.append("</svg>")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")


def gather_nwb_paths(inputs: Iterable[str]) -> list[Path]:
    """Expand user CLI inputs into a list of NWB file paths."""
    paths: list[Path] = []
    for item in inputs:
        entry = Path(item).expanduser().resolve()
        if entry.is_dir():
            paths.extend(sorted(entry.glob("*.nwb")))
        elif entry.suffix.lower() == ".nwb":
            paths.append(entry)
        else:
            raise ValueError(
                f"Unsupported input {entry}. Provide .nwb files or directories."
            )
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Detect saturated regions in NWB acquisitions and visualize them."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Paths to NWB files or directories containing NWB files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/saturation_plots"),
        help="Directory where SVG plots will be written.",
    )
    parser.add_argument(
        "--amp-factor",
        type=float,
        default=0.2,
        help="Absolute voltage threshold in volts (default 0.2 V = 200 mV). "
        "Values exceeding ±amp_factor are considered saturated.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=6000,
        help="Maximum number of plotted points per series (downsampling for SVG size).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    nwb_paths = gather_nwb_paths(args.inputs)
    if not nwb_paths:
        raise SystemExit("No NWB files found for the provided inputs.")

    for nwb_path in nwb_paths:
        print(f"Processing {nwb_path} ...")
        results = analyze_file(
            nwb_path=nwb_path,
            amp_factor=args.amp_factor,
            max_points=args.max_points,
        )

        # Create subject-specific output directory
        subject_dir = output_dir / nwb_path.stem
        subject_dir.mkdir(parents=True, exist_ok=True)

        # Separate clean and corrupted sweeps
        corrupted_results = [res for res in results if len(res.rectangles) > 0]
        clean_results = [res for res in results if len(res.rectangles) == 0]

        total_rects = sum(len(res.rectangles) for res in results)
        print(f"  • Sweeps processed: {len(results)}")
        print(f"  • Saturated segments detected: {total_rects}")
        print(f"  • Corrupted sweeps: {len(corrupted_results)}")
        print(f"  • Clean sweeps: {len(clean_results)}")

        # Generate one PNG per sweep (both clean and corrupted)
        for res in results:
            # Clean up series name for filename
            safe_name = res.series_name.replace("/", "_").replace(" ", "_")

            # Add suffix based on corruption status
            num_corrupted_regions = len(res.rectangles)
            if num_corrupted_regions == 0:
                suffix = "_clean"
            else:
                suffix = f"_corrupted_{num_corrupted_regions}"

            output_path = subject_dir / f"{safe_name}{suffix}.png"
            title = f"{res.series_name} ({suffix[1:]})"
            render_png([res], output_path=output_path, title=title)

        print(f"  • {len(results)} plots written to {subject_dir}/")


if __name__ == "__main__":
    main()
