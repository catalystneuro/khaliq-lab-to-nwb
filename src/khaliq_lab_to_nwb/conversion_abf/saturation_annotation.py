"""
Utilities for detecting and annotating voltage saturation in intracellular recordings.

This module provides functions to automatically detect saturated voltage segments
in NWB acquisition time series and annotate them in the invalid_times table.
"""

import numpy as np
from pynwb import NWBFile

from khaliq_lab_to_nwb.conversion_abf.detect_saturation import (
    detect_saturation_segments,
)


def add_saturation_annotations(
    nwbfile: NWBFile,
    voltage_threshold: float = 0.2,
) -> int:
    """
    Detect and annotate voltage saturation in all acquisition time series.

    This function scans all CurrentClampSeries in the NWB file's acquisition group,
    detects saturated voltage segments (exceeding ±voltage_threshold), and adds
    annotations to the invalid_times table with references to the affected time series.

    Only data that exists in the NWB file is annotated. Files that failed to load
    or parse are not included in annotations.

    Parameters
    ----------
    nwbfile : NWBFile
        NWB file with acquisition time series to scan for saturation.
    voltage_threshold : float, optional
        Absolute voltage threshold in volts (default 0.2 V = 200 mV).
        Voltages exceeding ±voltage_threshold are considered saturated.

    Returns
    -------
    int
        Number of saturation intervals detected and annotated.

    Notes
    -----
    The function adds custom columns to the invalid_times table:
    - time_series: Reference to the TimeSeries where saturation was detected
    - time_series_name: Name of the TimeSeries
    - voltage_min: Minimum voltage (V) during the saturated interval
    - voltage_max: Maximum voltage (V) during the saturated interval
    - reason: Set to "voltage_saturation_detected"
    - severity: Set to "high"

    Examples
    --------
    >>> from pynwb import NWBHDF5IO
    >>> with NWBHDF5IO('session.nwb', 'r+') as io:
    ...     nwbfile = io.read()
    ...     num_saturated = add_saturation_annotations(nwbfile)
    ...     print(f"Detected {num_saturated} saturation intervals")
    """
    # Detect saturation in all acquisition time series
    saturation_intervals = []

    for series_name, series in nwbfile.acquisition.items():
        # Get data and timing information
        data = np.array(series.data[:], dtype=np.float64)
        if series.timestamps is not None:
            start_time = float(series.timestamps[0])
        else:
            start_time = float(series.starting_time)
        rate = float(series.rate)

        # Detect saturation segments
        # Note: only voltage_threshold is actually used by detect_saturation_segments
        segments = detect_saturation_segments(
            data=data,
            rate=rate,
            start_time=start_time,
            amp_factor=voltage_threshold,
        )

        # Store segments with reference to the time series
        for start, stop, vmin, vmax in segments:
            saturation_intervals.append(
                {
                    "start_time": start,
                    "stop_time": stop,
                    "time_series": series,
                    "time_series_name": series_name,
                    "voltage_min": vmin,
                    "voltage_max": vmax,
                }
            )

    # Add invalid_times table if saturation detected
    # Only annotate data that actually exists in the NWB file
    if saturation_intervals:
        # Add custom columns to invalid_times table
        nwbfile.add_invalid_times_column(
            name="reason",
            description="Type of corruption or artifact detected in the data",
        )
        nwbfile.add_invalid_times_column(
            name="severity",
            description="Severity level of the corruption",
        )
        nwbfile.add_invalid_times_column(
            name="time_series",
            description="Reference to the TimeSeries where this invalid interval was detected",
        )
        nwbfile.add_invalid_times_column(
            name="time_series_name",
            description="Name of the TimeSeries where this invalid interval was detected",
        )
        nwbfile.add_invalid_times_column(
            name="voltage_min",
            description="Minimum voltage (V) during this saturated interval",
        )
        nwbfile.add_invalid_times_column(
            name="voltage_max",
            description="Maximum voltage (V) during this saturated interval",
        )

        # Add one entry per saturation interval
        for interval in saturation_intervals:
            nwbfile.add_invalid_time_interval(
                start_time=interval["start_time"],
                stop_time=interval["stop_time"],
                reason="voltage_saturation_detected",
                severity="high",
                time_series=interval["time_series"],
                time_series_name=interval["time_series_name"],
                voltage_min=interval["voltage_min"],
                voltage_max=interval["voltage_max"],
            )

    return len(saturation_intervals)
