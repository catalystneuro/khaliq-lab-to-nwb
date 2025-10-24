"""
Utilities for creating epoch annotations based on recording protocols.

This module provides functions to add protocol-based epochs to the NWB epochs table,
organizing recordings by their temporal structure (SF, CS, ADP protocol sequences).
"""

from pynwb import NWBFile


def add_protocol_epochs(
    nwbfile: NWBFile,
    files_with_recordings: list[dict],
) -> None:
    """
    Add protocol-based epochs to the NWB epochs table.

    Creates epochs that group recordings by protocol type (SF, CS, ADP) and cell,
    providing temporal organization of the experimental sequence. Each epoch represents
    one protocol recording session for a specific cell.

    Parameters
    ----------
    nwbfile : NWBFile
        NWB file to add epochs to.
    files_with_recordings : list[dict]
        List of file info dictionaries, each containing:
        - 'cell_id': cell identifier (e.g., "C1")
        - 'protocol': protocol abbreviation (SF, CS, ADP)
        - 'ais_status': AIS status string
        - 'recording_indices': list of intracellular recording indices
        - 'start_time': starting time of first sweep in this protocol
        - 'stop_time': stopping time of last sweep in this protocol

    Notes
    -----
    The epochs table will have custom columns added:
    - protocol: Protocol abbreviation (SF, CS, ADP)
    - cell_id: Cell identifier
    - ais_status: Axon initial segment component status
    - tags: List of descriptive tags for the protocol type

    Protocol tags:
    - SF: ["spontaneous_firing", "baseline"]
    - CS: ["current_steps", "excitability_test"]
    - ADP: ["after_depolarization", "post_stimulus"]

    Examples
    --------
    >>> # After adding all recordings to nwbfile
    >>> add_protocol_epochs(nwbfile, files_with_recordings)
    >>> # Now can query epochs
    >>> sf_epochs = [e for e in nwbfile.epochs.to_dataframe().iterrows()
    ...              if e[1]['protocol'] == 'SF']
    """
    # Protocol-specific tags for semantic annotation
    protocol_tags = {
        "SF": ["spontaneous_firing", "baseline"],
        "CS": ["current_steps", "excitability_test"],
        "ADP": ["after_depolarization", "post_stimulus"],
    }

    # Protocol descriptions
    protocol_descriptions = {
        "SF": "Spontaneous Firing: Baseline intrinsic firing pattern without current injection",
        "CS": "Current Steps: Response to stepwise current injection of varying amplitudes",
        "ADP": "After Depolarization: Post-stimulus depolarization following current injection",
    }

    # Add custom columns to epochs table
    nwbfile.add_epoch_column(
        name="protocol",
        description="Protocol abbreviation (SF=Spontaneous Firing, CS=Current Steps, ADP=After Depolarization)",
    )
    nwbfile.add_epoch_column(
        name="cell_id",
        description="Cell identifier (e.g., C1, C2, ...)",
    )
    nwbfile.add_epoch_column(
        name="ais_status",
        description="Axon initial segment component status",
    )

    # Add one epoch per protocol recording
    for file_info in files_with_recordings:
        cell_id = file_info["cell_id"]
        protocol = file_info["protocol"]
        ais_status = file_info["ais_status"]
        start_time = file_info["start_time"]
        stop_time = file_info["stop_time"]

        # Create epoch
        nwbfile.add_epoch(
            start_time=start_time,
            stop_time=stop_time,
            tags=protocol_tags.get(protocol, [protocol]),
            protocol=protocol,
            cell_id=cell_id,
            ais_status=ais_status,
        )
