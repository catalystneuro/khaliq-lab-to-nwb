"""
Functions for adding intracellular current-clamp electrophysiology data to NWB files using Neo.

This module provides direct control over how voltage (response) and current (stimulus)
signals are mapped from ABF files to NWB's intracellular ephys structures.
"""

from neo.rawio import AxonRawIO
from pynwb import NWBFile
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries


def add_icephys_data_from_neo_reader(
    nwbfile: NWBFile,
    neo_reader: AxonRawIO,
    electrode_name: str = "electrode0",
) -> NWBFile:
    """
    Add intracellular current-clamp electrophysiology data from a Neo AxonRawIO reader to an NWBFile.

    This function reads ABF data using a Neo reader instance and properly maps the signals based on units:
    - Voltage units (mV, V) → CurrentClampSeries (response - what is MEASURED)
    - Current units (pA, nA, A) → CurrentClampStimulusSeries (stimulus - what is CONTROLLED)

    Each sweep in the ABF file becomes a separate stimulus-response pair in NWB.

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file object to add data to. Must already have an intracellular
        electrode defined with name matching `electrode_name`.
    neo_reader : AxonRawIO
        Neo AxonRawIO reader instance with header already parsed (parse_header() called).
    electrode_name : str, default: "electrode0"
        Name of the electrode in the NWBFile to associate with these recordings.
        This electrode must already exist in nwbfile.icephys_electrodes.

    Returns
    -------
    NWBFile
        The NWBFile with intracellular data added.

    Raises
    ------
    KeyError
        If the specified electrode_name is not found in the NWBFile.
    ValueError
        If no voltage channel is found.

    Notes
    -----
    **Channel Identification (by units only):**
    - Voltage channel: mV or V → response
    - Current channel: pA, nA, or A → stimulus
    - Other units (C, etc.) → ignored

    **Unit Conversions to SI base units:**
    - Voltage: mV → V (divide by 1000)
    - Current: pA → A (divide by 1e12), nA → A (divide by 1e9)

    **Protocol Handling:**
    - SF (Spontaneous Firing): Only voltage channel, no current
    - CS/ADP protocols: Both voltage and current channels

    Examples
    --------
    >>> from pathlib import Path
    >>> from neo.rawio import AxonRawIO
    >>> from pynwb import NWBFile
    >>> from datetime import datetime
    >>> from zoneinfo import ZoneInfo
    >>>
    >>> # Initialize Neo reader
    >>> abf_path = Path("recording.abf")
    >>> reader = AxonRawIO(filename=abf_path)
    >>> reader.parse_header()
    >>>
    >>> # Get session start time from reader
    >>> session_start_time = reader.raw_annotations['blocks'][0]['rec_datetime']
    >>>
    >>> # Create NWB file with electrode
    >>> nwbfile = NWBFile(
    ...     session_description="Test session",
    ...     identifier="test_001",
    ...     session_start_time=session_start_time
    ... )
    >>> device = nwbfile.create_device(name="MultiClamp700B")
    >>> electrode = nwbfile.create_icephys_electrode(
    ...     name="electrode0",
    ...     description="Patch pipette",
    ...     device=device
    ... )
    >>>
    >>> # Add current-clamp data from Neo reader
    >>> add_icephys_data_from_neo_reader(
    ...     nwbfile=nwbfile,
    ...     neo_reader=reader,
    ...     electrode_name="electrode0"
    ... )
    """
    # Get the electrode from the NWBFile
    if electrode_name not in nwbfile.icephys_electrodes:
        raise KeyError(
            f"Electrode '{electrode_name}' not found in NWBFile. "
            f"Available electrodes: {list(nwbfile.icephys_electrodes.keys())}"
        )
    electrode = nwbfile.icephys_electrodes[electrode_name]

    # Identify voltage and current channels by units
    signal_channels = neo_reader.header["signal_channels"]
    voltage_channel_index = None
    current_channel_index = None

    for index, channel in enumerate(signal_channels):
        units = channel["units"]

        # Voltage channel: mV or V
        if units in ("mV", "V"):
            voltage_channel_index = index

        # Current channel: pA, nA, or A
        elif units in ("pA", "nA", "A"):
            current_channel_index = index

    # Voltage channel must be present
    if voltage_channel_index is None:
        raise ValueError(
            f"Could not identify voltage channel (expected units: mV or V). "
            f"Found channels: {[(ch['name'], ch['units']) for ch in signal_channels]}"
        )

    # Get channel info
    voltage_channel = signal_channels[voltage_channel_index]
    sampling_rate = float(voltage_channel["sampling_rate"])

    # Check if current channel is available
    has_current_channel = current_channel_index is not None
    if has_current_channel:
        current_channel = signal_channels[current_channel_index]

    # Process each sweep (segment)
    nb_segments = neo_reader.header["nb_segment"][0]

    for segment_index in range(nb_segments):
        # Get timing information
        t_start = neo_reader._segment_t_start(block_index=0, seg_index=segment_index)

        # Read voltage data
        voltage_raw = neo_reader.get_analogsignal_chunk(
            block_index=0,
            seg_index=segment_index,
            i_start=None,
            i_stop=None,
            channel_indexes=[voltage_channel_index],
        )
        voltage_data = (
            voltage_raw.astype("float32") * voltage_channel["gain"]
            + voltage_channel["offset"]
        ).squeeze()

        # Convert voltage to volts
        if voltage_channel["units"] == "mV":
            voltage_data = voltage_data / 1000.0

        # Read current data if available
        if has_current_channel:
            current_raw = neo_reader.get_analogsignal_chunk(
                block_index=0,
                seg_index=segment_index,
                i_start=None,
                i_stop=None,
                channel_indexes=[current_channel_index],
            )
            current_data = (
                current_raw.astype("float32") * current_channel["gain"]
                + current_channel["offset"]
            ).squeeze()

            # Convert current to amperes
            if current_channel["units"] == "pA":
                current_data = current_data / 1e12
            elif current_channel["units"] == "nA":
                current_data = current_data / 1e9

        # Create unique names for this sweep
        sweep_number = segment_index
        response_name = f"CurrentClampSeries_{sweep_number:03d}"

        # Add stimulus if current channel is present
        if has_current_channel:
            stimulus_name = f"CurrentClampStimulusSeries_{sweep_number:03d}"

            stimulus_series = CurrentClampStimulusSeries(
                name=stimulus_name,
                data=current_data,
                starting_time=t_start,
                rate=sampling_rate,
                electrode=electrode,
                gain=1.0,
                unit="amperes",
                description=f"Current stimulus for sweep {sweep_number}",
            )
            nwbfile.add_stimulus(stimulus_series)
        else:
            # No current injection (spontaneous firing)
            stimulus_series = None

        # Add voltage response
        response_series = CurrentClampSeries(
            name=response_name,
            data=voltage_data,
            starting_time=t_start,
            rate=sampling_rate,
            electrode=electrode,
            gain=1.0,
            unit="volts",
            stimulus_description=stimulus_series.name if stimulus_series else "none",
            description=f"Voltage response for sweep {sweep_number}",
        )
        nwbfile.add_acquisition(response_series)

        # Add intracellular recording to link stimulus and response
        nwbfile.add_intracellular_recording(
            electrode=electrode,
            stimulus=stimulus_series,
            response=response_series,
        )

    return nwbfile
