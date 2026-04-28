from typing import Literal

import numpy as np
from neo.rawio import AxonRawIO
from pynwb import NWBFile
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries


def _identify_signal_channels(
    neo_reader: AxonRawIO,
) -> tuple[int, int | None]:
    """Locate the voltage and current channels in a Neo AxonRawIO reader by units."""
    signal_channels = neo_reader.header["signal_channels"]
    voltage_channel_index = None
    current_channel_index = None

    for index, channel in enumerate(signal_channels):
        units = channel["units"]
        if units in ("mV", "V"):
            voltage_channel_index = index
        elif units in ("pA", "nA", "A"):
            current_channel_index = index

    if voltage_channel_index is None:
        raise ValueError(
            f"Could not identify voltage channel (expected units: mV or V). "
            f"Found channels: {[(ch['name'], ch['units']) for ch in signal_channels]}"
        )

    return voltage_channel_index, current_channel_index


def _read_segment(
    neo_reader: AxonRawIO,
    segment_index: int,
    channel_index: int,
    target_units: Literal["volts", "amperes"],
) -> np.ndarray:
    """Read one segment from a channel and convert to SI units."""
    channel = neo_reader.header["signal_channels"][channel_index]
    raw = neo_reader.get_analogsignal_chunk(
        block_index=0,
        seg_index=segment_index,
        i_start=None,
        i_stop=None,
        channel_indexes=[channel_index],
    )
    data = (raw.astype("float32") * channel["gain"] + channel["offset"]).squeeze()

    units = channel["units"]
    if target_units == "volts":
        if units == "mV":
            data = data / 1000.0
    elif target_units == "amperes":
        if units == "pA":
            data = data / 1e12
        elif units == "nA":
            data = data / 1e9

    return data


def add_icephys_data_from_neo_readers(
    nwbfile: NWBFile,
    neo_readers: list[AxonRawIO],
    time_offsets: list[float],
    electrode_name: str,
    cell_id: str,
    ais_string: Literal["AIS", "NOAIS"],
    protocol: Literal["SF", "CS", "ADP"],
) -> tuple[list[int], float, float]:
    """
    Add a single concatenated current-clamp recording for one (cell, protocol) group.

    All sweeps from all provided ABF files are concatenated into one CurrentClampSeries
    (response) and, when a current channel is present, one CurrentClampStimulusSeries
    (stimulus). A timestamps array preserves the real wall-clock spacing between sweeps
    and between files. Each individual sweep is exposed as a region of these series via
    rows in the IntracellularRecordingsTable, using the start_index / index_count
    parameters of add_intracellular_recording.

    The two series share timestamps via a soft link: pynwb stores the timestamps array
    on the response series and the stimulus series references it.

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add data to. Must contain an electrode named ``electrode_name``.
    neo_readers : list[AxonRawIO]
        Neo readers for all ABF files belonging to this (cell, protocol) group.
        Header must already be parsed.
    time_offsets : list[float]
        Per-reader offset (seconds) from the session start time. Same length as ``neo_readers``.
    electrode_name : str
        Name of the icephys electrode to attach the recordings to.
    cell_id : str
        Cell identifier (e.g. "C1"). Used to build series names.
    ais_string : {"AIS", "NOAIS"}
        AIS status string used in series naming.
    protocol : {"SF", "CS", "ADP"}
        Protocol abbreviation used in series naming.

    Returns
    -------
    recording_indices : list[int]
        IntracellularRecordingsTable row indices, one per sweep, in the order sweeps
        were added (file order, then segment order within each file).
    start_time : float
        Absolute start time (s) of the first sample of the first sweep, including offset.
    stop_time : float
        Absolute stop time (s) of the last sample of the last sweep, including offset.
    """
    if len(neo_readers) != len(time_offsets):
        raise ValueError("neo_readers and time_offsets must have the same length")
    if not neo_readers:
        raise ValueError("neo_readers must contain at least one reader")

    if electrode_name not in nwbfile.icephys_electrodes:
        raise KeyError(
            f"Electrode '{electrode_name}' not found in NWBFile. "
            f"Available electrodes: {list(nwbfile.icephys_electrodes.keys())}"
        )
    electrode = nwbfile.icephys_electrodes[electrode_name]

    voltage_channel_index, current_channel_index = _identify_signal_channels(
        neo_readers[0]
    )
    has_current_channel = current_channel_index is not None

    voltage_arrays: list[np.ndarray] = []
    current_arrays: list[np.ndarray] = []
    timestamp_arrays: list[np.ndarray] = []
    sweep_regions: list[tuple[int, int]] = []

    cumulative_index = 0
    for reader, time_offset in zip(neo_readers, time_offsets):
        v_idx, c_idx = _identify_signal_channels(reader)
        if (c_idx is None) != (not has_current_channel):
            raise ValueError(
                "Inconsistent channel layout across files in the same (cell, protocol) "
                "group: some files have a current channel and others do not."
            )

        sampling_rate = float(reader.header["signal_channels"][v_idx]["sampling_rate"])
        nb_segments = reader.header["nb_segment"][0]

        for segment_index in range(nb_segments):
            t_start = reader._segment_t_start(block_index=0, seg_index=segment_index)

            voltage_data = _read_segment(reader, segment_index, v_idx, "volts")
            n_samples = voltage_data.size

            timestamps = (t_start + time_offset) + np.arange(
                n_samples, dtype=np.float64
            ) / sampling_rate

            voltage_arrays.append(voltage_data)
            timestamp_arrays.append(timestamps)

            if has_current_channel:
                current_data = _read_segment(reader, segment_index, c_idx, "amperes")
                current_arrays.append(current_data)

            sweep_regions.append((cumulative_index, n_samples))
            cumulative_index += n_samples

    all_voltage = np.concatenate(voltage_arrays)
    all_timestamps = np.concatenate(timestamp_arrays)
    all_current = np.concatenate(current_arrays) if has_current_channel else None

    cell_number_str = cell_id[1:] if cell_id else "00"
    ais_part = "Ais" if ais_string == "AIS" else "Noais"
    response_name = (
        f"CurrentClampSeries{ais_part}Protocol{protocol}Cell{cell_number_str:0>2}"
    )
    stimulus_name = f"CurrentClampStimulusSeries{ais_part}Protocol{protocol}Cell{cell_number_str:0>2}"

    response_series = CurrentClampSeries(
        name=response_name,
        data=all_voltage,
        timestamps=all_timestamps,
        electrode=electrode,
        gain=1.0,
        unit="volts",
        stimulus_description=stimulus_name if has_current_channel else "none",
        description=(
            f"Concatenated voltage response across all sweeps of the {protocol} "
            f"protocol for cell {cell_id}. Individual sweeps are exposed as rows of "
            f"the IntracellularRecordingsTable, each pointing to a region of this "
            f"series via response_start_index and response_index_count."
        ),
    )
    nwbfile.add_acquisition(response_series)

    stimulus_series = None
    if has_current_channel:
        stimulus_series = CurrentClampStimulusSeries(
            name=stimulus_name,
            data=all_current,
            timestamps=response_series,
            electrode=electrode,
            gain=1.0,
            unit="amperes",
            description=(
                f"Concatenated current stimulus across all sweeps of the {protocol} "
                f"protocol for cell {cell_id}. Individual sweeps are exposed as rows "
                f"of the IntracellularRecordingsTable, each pointing to a region of "
                f"this series via stimulus_start_index and stimulus_index_count."
            ),
        )
        nwbfile.add_stimulus(stimulus_series)

    recording_indices: list[int] = []
    for start_idx, count in sweep_regions:
        if stimulus_series is not None:
            ir_index = nwbfile.add_intracellular_recording(
                electrode=electrode,
                stimulus=stimulus_series,
                stimulus_start_index=start_idx,
                stimulus_index_count=count,
                response=response_series,
                response_start_index=start_idx,
                response_index_count=count,
            )
        else:
            ir_index = nwbfile.add_intracellular_recording(
                electrode=electrode,
                response=response_series,
                response_start_index=start_idx,
                response_index_count=count,
            )
        recording_indices.append(ir_index)

    start_time = float(all_timestamps[0])
    stop_time = float(all_timestamps[-1])
    return recording_indices, start_time, stop_time
