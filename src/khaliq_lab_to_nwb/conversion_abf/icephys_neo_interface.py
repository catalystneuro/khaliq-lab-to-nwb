from typing import Literal

from neo.rawio import AxonRawIO
from pynwb import NWBFile
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries


def add_icephys_data_from_neo_reader(
    nwbfile: NWBFile,
    neo_reader: AxonRawIO,
    electrode_name: str = "electrode0",
    time_offset: float = 0.0,
    cell_id: str = "",
    sweep_offset: int = 0,
    ais_string: Literal["AIS", "NOAIS"] = "AIS",
    protocol: Literal["SF", "CS", "ADP"] = "SF",
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
    time_offset : float, default: 0.0
        Time offset in seconds to add to all timestamps. Used when combining multiple
        recordings into a single NWB file to adjust timestamps relative to a common
        reference time (typically the earliest recording in the session).
    cell_id : str, default: ""
        Cell identifier to include in sweep names. Used when multiple cells are
        recorded in the same NWB file to distinguish sweeps from different cells.
    sweep_offset : int, default: 0
        Offset to add to sweep numbers. Used when multiple files for the same cell
        are combined to avoid naming conflicts. Each file's sweeps will be numbered
        starting from this offset.
    ais_string : str, default: ""
        AIS status string ("AIS" or "NOAIS") for naming series.
    protocol : str, default: ""
        Protocol abbreviation (SF, CS, ADP) for naming series.

    Returns
    -------
    list[int]
        List of intracellular recording indices (one per sweep added).

    Raises
    ------
    KeyError
        If the specified electrode_name is not found in the NWBFile.
    ValueError
        If no voltage channel is found.
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
    recording_indices = []

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
        # Apply sweep offset to avoid conflicts when combining multiple files
        sweep_number = segment_index + sweep_offset

        # Extract cell number from cell_id (e.g., "C14" -> "14")
        cell_number_str = cell_id[1:] if cell_id else "00"

        # Format: CurrentClampSeriesAisProtocolSfCell01Sweep022
        # Components: {Ais|Noais}Protocol{XX}Cell{XX}Sweep{XXX}
        ais_part = "Ais" if ais_string == "AIS" else "Noais"
        response_name = f"CurrentClampSeries{ais_part}Protocol{protocol}Cell{cell_number_str:0>2}Sweep{sweep_number:03d}"

        # Adjust timestamp relative to reference time
        adjusted_t_start = t_start + time_offset

        # Add stimulus if current channel is present
        if has_current_channel:
            stimulus_name = f"CurrentClampStimulusSeries{ais_part}Protocol{protocol}Cell{cell_number_str:0>2}Sweep{sweep_number:03d}"

            stimulus_series = CurrentClampStimulusSeries(
                name=stimulus_name,
                data=current_data,
                starting_time=adjusted_t_start,
                rate=sampling_rate,
                electrode=electrode,
                gain=1.0,
                unit="amperes",
                description=f"Current stimulus for sweep {sweep_number}"
                + (f" (cell {cell_id})" if cell_id else ""),
            )
            nwbfile.add_stimulus(stimulus_series)
        else:
            # No current injection (spontaneous firing)
            stimulus_series = None

        # Add voltage response
        response_series = CurrentClampSeries(
            name=response_name,
            data=voltage_data,
            starting_time=adjusted_t_start,
            rate=sampling_rate,
            electrode=electrode,
            gain=1.0,
            unit="volts",
            stimulus_description=stimulus_series.name if stimulus_series else "none",
            description=f"Voltage response for sweep {sweep_number}"
            + (f" (cell {cell_id})" if cell_id else ""),
        )
        nwbfile.add_acquisition(response_series)

        # Add intracellular recording to link stimulus and response
        ir_index = nwbfile.add_intracellular_recording(
            electrode=electrode,
            stimulus=stimulus_series,
            response=response_series,
        )
        recording_indices.append(ir_index)

    return recording_indices
