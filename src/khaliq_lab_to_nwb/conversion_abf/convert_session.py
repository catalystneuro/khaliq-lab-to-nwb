from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd
from neo.rawio import AxonRawIO
from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile
from pynwb.file import Subject
from tqdm import tqdm

from khaliq_lab_to_nwb.conversion_abf.icephys_neo_interface import (
    add_icephys_data_from_neo_reader,
)


def check_corruption_status(file_path: Path) -> bool:
    """
    Check if an ABF file is marked as corrupted in the corruption status CSV.

    This function looks up the given file path in a CSV database that tracks
    which ABF files have been identified as corrupted or containing artifacts.
    Files not found in the CSV are assumed to be clean.

    Parameters
    ----------
    file_path : Path
        Path to the ABF file to check

    Returns
    -------
    bool
        True if file is marked as corrupted, False if clean or not found in CSV

    Notes
    -----
    - The corruption CSV file is located at assets/file_corruption_status.csv
    - CSV contains columns: "file_path" (absolute path) and "is_corrupted" (boolean)
    - Paths are normalized and resolved for reliable comparison
    - Unknown files (not in CSV) are assumed to be clean
    - This is used to add invalid_times annotations to NWB files
    """
    # Get the corruption CSV path relative to this script
    corruption_csv_path = (
        Path(__file__).parent / "assets" / "file_corruption_status.csv"
    )
    df = pd.read_csv(corruption_csv_path)

    # Normalize paths for comparison
    file_path_resolved = file_path.resolve()

    for _, row in df.iterrows():
        csv_path = Path(row["file_path"]).resolve()
        if csv_path == file_path_resolved:
            return row["is_corrupted"]

    # If not found in CSV, assume not corrupted
    return False


def extract_cell_number_from_path(file_path: Path) -> int:
    """
    Extract cell number from the ABF file path structure.

    The Khaliq Lab data follows a directory naming convention where each cell
    is stored in a folder named "C#" (e.g., "C1", "C2", "C3"). This function
    parses the path to extract that number.

    Parameters
    ----------
    file_path : Path
        Path to the ABF file, expected to contain a directory component
        matching the pattern "C#" where # is one or more digits

    Returns
    -------
    int
        The cell number extracted from the path

    Raises
    ------
    ValueError
        If no directory matching the "C#" pattern is found in the path
    """
    # The cell number is typically in a directory like "C1", "C2", etc.
    parts = file_path.parts
    for part in parts:
        if part.startswith("C") and part[1:].isdigit():
            return int(part[1:])

    raise ValueError(f"Could not extract cell number from path: {file_path}")


def get_cells_for_subject(subject_id: str, session_metadata_path: Path) -> list[dict]:
    """
    Get all cell information for a specific subject from metadata CSV.

    Parameters
    ----------
    subject_id : str
        Subject identifier (e.g., "#1", "#2")
    session_metadata_path : Path
        Path to the session metadata CSV file

    Returns
    -------
    list[dict]
        List of dictionaries, each containing metadata for one cell belonging to the subject
    """
    df = pd.read_csv(session_metadata_path)
    subject_cells = df[df["Animal ID"] == subject_id]

    if subject_cells.empty:
        raise ValueError(f"No cells found for subject {subject_id}")

    cells = []
    for _, row in subject_cells.iterrows():
        age_str = row["Animal Age (Years)"]
        age_years = None if age_str == "N.A." else float(age_str)

        cells.append(
            {
                "date": row["Date"],
                "cell_number": int(row["Cell #"]),
                "neuron_type": row["Neuron Type"],
                "anatomical_region": row["Anatomical Region"],
                "animal_id": row["Animal ID"],
                "animal_species": row["Animal Species"],
                "sex": row["Sex"],
                "age_years": age_years,
            }
        )

    return cells


def collect_subject_files(
    subject_id: str, session_metadata_path: Path, data_folder_path: Path
) -> list[dict]:
    """
    Collect all ABF files for a specific subject.

    Parameters
    ----------
    subject_id : str
        Subject identifier (e.g., "#1", "#2")
    session_metadata_path : Path
        Path to the session metadata CSV file
    data_folder_path : Path
        Root folder containing electrophysiology recordings

    Returns
    -------
    list[dict]
        List of dictionaries, each containing:
        - "file_path": Path to ABF file
        - "cell_number": Cell number
        - "cell_id": Cell ID string (e.g., "C1")
        - "ais_status": "with AIS component" or "without AIS component"
        - "protocol": Protocol type (SF, CS, ADP)
        - "rec_datetime": Recording datetime from ABF file
    """
    # Get all cells for this subject
    cells = get_cells_for_subject(subject_id, session_metadata_path)
    cell_numbers = [cell["cell_number"] for cell in cells]

    # Search for ABF files
    all_files = []
    for ais_folder in ["Cells with AIS component", "Cells without AIS component"]:
        ais_path = data_folder_path / ais_folder
        if not ais_path.exists():
            continue

        ais_status = (
            "with AIS component"
            if ais_folder == "Cells with AIS component"
            else "without AIS component"
        )

        # Find all ABF files
        for abf_file in ais_path.glob("**/*.abf"):
            # Extract cell number from path
            cell_number = extract_cell_number_from_path(abf_file)

            # Skip if not belonging to this subject
            if cell_number not in cell_numbers:
                continue

            # Skip the misplaced C3 file
            if cell_number == 3 and abf_file.name == "21610007.abf":
                continue

            # Get protocol from parent directory
            protocol = abf_file.parent.name

            all_files.append(
                {
                    "file_path": abf_file,
                    "cell_number": cell_number,
                    "cell_id": f"C{cell_number}",
                    "ais_status": ais_status,
                    "protocol": protocol,
                }
            )

    return all_files


# ============================================================================
# MAIN CONVERSION WORKFLOW
# ============================================================================


def convert_subject(
    subject_id: str,
    session_metadata_path: Path | str,
    data_folder_path: Path | str,
    output_folder_path: Path | str,
) -> Path:
    """
    Convert all recordings for a single subject to one NWB file.

    This function collects all ABF files for a subject (across all cells and protocols),
    creates a single NWB file with multiple electrodes, and adjusts all timestamps
    relative to the earliest recording.

    Parameters
    ----------
    subject_id : str
        Subject identifier (e.g., "#1", "#2")
    session_metadata_path : Path or str
        Path to the session metadata CSV file
    data_folder_path : Path or str
        Root folder containing electrophysiology recordings
    output_folder_path : Path or str
        Path to folder where NWB file will be saved

    Returns
    -------
    Path
        Path to the created NWB file
    """
    session_metadata_path = Path(session_metadata_path)
    data_folder_path = Path(data_folder_path)
    output_folder_path = Path(output_folder_path)

    # Collect all files for this subject
    all_files = collect_subject_files(
        subject_id, session_metadata_path, data_folder_path
    )

    if not all_files:
        raise ValueError(f"No ABF files found for subject {subject_id}")

    # Open all files once to get timing information
    # Skip files that cannot be parsed (corrupted header)
    files_with_timing = []
    corrupted_files = []
    unparsable_files = []

    for file_info in all_files:
        file_path = file_info["file_path"]
        is_corrupted = check_corruption_status(file_path)

        try:
            neo_reader = AxonRawIO(filename=str(file_path))
            neo_reader.parse_header()
            rec_datetime = neo_reader.raw_annotations["blocks"][0]["rec_datetime"]

            files_with_timing.append(
                {
                    **file_info,
                    "rec_datetime": rec_datetime,
                    "neo_reader": neo_reader,
                }
            )

            if is_corrupted:
                corrupted_files.append((file_path.name, file_info["cell_id"]))
        except Exception as e:
            unparsable_files.append((file_path.name, file_info["cell_id"]))

    # Sort files by recording time to ensure chronological processing
    files_with_timing.sort(key=lambda f: f["rec_datetime"])

    # Find earliest recording time
    earliest_time = files_with_timing[0]["rec_datetime"]
    all_files = files_with_timing  # Use the updated list

    # Add timezone info (NIH Bethesda, MD is in US Eastern Time)
    session_start_time = earliest_time.replace(tzinfo=ZoneInfo("America/New_York"))

    # Load subject metadata from first cell
    cells_metadata = get_cells_for_subject(subject_id, session_metadata_path)
    first_cell = cells_metadata[0]

    # Convert species
    species_mapping = {
        "Rhesus Macaque": "Macaca mulatta",
        "Rhesus macaque": "Macaca mulatta",
    }
    species_latin = species_mapping.get(
        first_cell["animal_species"], first_cell["animal_species"]
    )

    # Convert age from years to days in ISO 8601 format
    if first_cell["age_years"] is not None:
        age_iso = f"P{first_cell['age_years'] * 365:.0f}D"
    else:
        age_iso = None

    # Count unique cells and protocols
    unique_cells = sorted(set(f["cell_id"] for f in all_files))
    unique_protocols = sorted(set(f["protocol"] for f in all_files))
    ais_statuses = sorted(set(f["ais_status"] for f in all_files))

    # Create Subject
    subject = Subject(
        subject_id=subject_id,
        species=species_latin,
        sex=first_cell["sex"],
        age=age_iso,
        description=f"Rhesus macaque subject with intracellular recordings from {first_cell['neuron_type']} neurons in {first_cell['anatomical_region']}",
    )

    # Protocol description - same for all subjects
    protocol_narrative = (
        "The experimental protocol consisted of three sequential current-clamp recordings "
        "to characterize intrinsic firing properties of dopaminergic neurons. "
        "First, spontaneous firing (SF) was recorded to establish baseline intrinsic firing "
        "patterns without current injection. "
        "Second, current steps (CS) protocols tested neuronal response to stepwise current "
        "injection of varying amplitudes. "
        "Third, after-depolarization (ADP) measurements examined post-stimulus depolarization "
        "following current injection. "
        "Recordings were performed on cells both with and without axon initial segment (AIS) "
        "components to investigate the role of AIS morphology in intrinsic excitability."
    )

    nwbfile = NWBFile(
        session_description=(
            f"Multi-cell intracellular current-clamp recordings from {len(unique_cells)} "
            f"{first_cell['neuron_type']} neuron(s) in {first_cell['anatomical_region']} "
            f"of subject {subject_id}"
        ),
        identifier=f"subject_{subject_id.replace('#', '')}",
        session_start_time=session_start_time,
        experimenter=["Khaliq, Zayd", "Sansalone, Lorenze"],
        lab="Cellular Neurophysiology Section",
        institution=(
            "National Institute of Neurological Disorders and Stroke, National Institutes of Health"
        ),
        experiment_description=(
            f"Intracellular patch-clamp electrophysiology recordings from dopaminergic neurons "
            f"in the substantia nigra pars compacta to investigate intrinsic firing properties "
            f"and the role of axon initial segment (AIS) morphology in neuronal excitability. "
            f"{protocol_narrative}"
        ),
        protocol=protocol_narrative,
        keywords=[
            first_cell["neuron_type"],
            "substantia nigra pars compacta",
            "patch-clamp",
            "intracellular",
            "current-clamp",
            "axon initial segment",
        ]
        + unique_protocols
        + ais_statuses,
        subject=subject,
    )

    # Step 6: Create device (shared across all electrodes)
    device = nwbfile.create_device(
        name="MultiClamp700B",
        description="Molecular Devices MultiClamp 700B amplifier",
        manufacturer="Molecular Devices",
    )

    # Create electrodes for each unique cell
    cell_to_electrode = {}

    # Group files by cell to get AIS status for each
    from collections import defaultdict

    cell_info = defaultdict(set)
    for f in all_files:
        cell_info[f["cell_id"]].add(f["ais_status"])

    for cell_id in unique_cells:
        # Get AIS status for this cell (should be consistent across all files for this cell)
        ais_status_set = cell_info[cell_id]
        if len(ais_status_set) > 1:
            raise ValueError(
                f"Cell {cell_id} has inconsistent AIS status: {ais_status_set}"
            )
        ais_status = list(ais_status_set)[0]

        cell_number = int(cell_id[1:])  # Extract number from "C#"
        electrode_name = f"IntracellularElectrodeCell{cell_number:02d}"

        # Map to standard brain atlas terminology
        # Substantia Nigra -> substantia nigra pars compacta (SNc) for dopaminergic neurons
        location_mapping = {
            "Substantia Nigra": "substantia nigra pars compacta",
        }
        location_standard = location_mapping.get(
            first_cell["anatomical_region"], first_cell["anatomical_region"].lower()
        )

        nwbfile.create_icephys_electrode(
            name=electrode_name,
            description=(
                f"Intracellular electrode for {first_cell['neuron_type']} neuron "
                f"(cell {cell_id}, {ais_status})"
            ),
            device=device,
            location=location_standard,
            cell_id=f"{first_cell['neuron_type']}Cell{cell_number:02d}",
        )

        cell_to_electrode[cell_id] = electrode_name

    # Add data from all files
    # Track sweep counters for each cell to avoid naming conflicts
    # Track files that fail during data reading
    cell_sweep_counters = defaultdict(int)
    data_read_failures = []

    for file_info in all_files:
        file_path = file_info["file_path"]
        cell_id = file_info["cell_id"]
        electrode_name = cell_to_electrode[cell_id]

        # Calculate time offset relative to session start
        file_rec_time = file_info["rec_datetime"]
        time_offset = (file_rec_time - earliest_time).total_seconds()

        # Use the already-opened neo_reader
        neo_reader = file_info["neo_reader"]

        # Get the current sweep counter for this cell
        sweep_offset = cell_sweep_counters[cell_id]

        # Determine AIS status string for naming
        ais_string = (
            "AIS" if file_info["ais_status"] == "with AIS component" else "NOAIS"
        )
        protocol = file_info["protocol"]

        try:
            add_icephys_data_from_neo_reader(
                nwbfile=nwbfile,
                neo_reader=neo_reader,
                electrode_name=electrode_name,
                time_offset=time_offset,
                cell_id=cell_id,
                sweep_offset=sweep_offset,
                ais_string=ais_string,
                protocol=protocol,
            )

            # Update sweep counter for this cell
            nb_segments = neo_reader.header["nb_segment"][0]
            cell_sweep_counters[cell_id] += nb_segments
        except (ValueError, OSError) as e:
            # File data cannot be read (e.g., mmap errors)
            data_read_failures.append((file_path.name, cell_id))

    # Add corruption annotations if any corrupted, unparsable, or failed files
    if corrupted_files or unparsable_files or data_read_failures:
        nwbfile.add_invalid_times_column(
            name="reason",
            description="Type of corruption or artifact detected in the data",
        )
        nwbfile.add_invalid_times_column(
            name="severity", description="Severity level of the corruption"
        )

        # Add a placeholder interval for each corrupted file (timing issues)
        for file_name, cell_id in corrupted_files:
            nwbfile.add_invalid_time_interval(
                start_time=0.0,
                stop_time=0.0,
                reason=f"file_marked_as_corrupted: {file_name} (cell {cell_id})",
                severity="unknown",
            )

        # Add entries for unparsable files (corrupted header)
        for file_name, cell_id in unparsable_files:
            nwbfile.add_invalid_time_interval(
                start_time=0.0,
                stop_time=0.0,
                reason=f"file_header_corrupted_cannot_parse: {file_name} (cell {cell_id})",
                severity="critical",
            )

        # Add entries for files that failed during data reading
        for file_name, cell_id in data_read_failures:
            nwbfile.add_invalid_time_interval(
                start_time=0.0,
                stop_time=0.0,
                reason=f"file_data_read_failed: {file_name} (cell {cell_id})",
                severity="critical",
            )

    # Write NWB file
    output_folder_path.mkdir(parents=True, exist_ok=True)

    output_filename = f"subject_{subject_id.replace('#', '')}.nwb"
    nwbfile_path = output_folder_path / output_filename

    configure_and_write_nwbfile(
        nwbfile=nwbfile, nwbfile_path=nwbfile_path, backend="hdf5"
    )

    return nwbfile_path


def main():
    """
    Main entry point for subject-level ABF to NWB conversion.

    This function converts all recordings for a single subject into one NWB file.
    Select the subject by modifying the SUBJECT_INDEX variable below.

    Available subjects:
    - Index 0: Subject #1 (5 cells)
    - Index 1: Subject #2 (29 cells)
    - Index 2: Subject #3 (2 cells)
    - Index 3: Subject #4 (3 cells)
    - Index 4: Subject #5 (1 cell)
    - Index 5: Subject #6 (1 cell)
    """

    # ============================================================================
    # PATHS
    # ============================================================================

    # Base paths
    repo_root_path = Path(__file__).parent.parent.parent.parent
    data_folder_path = repo_root_path / "data" / "Electrophysiology recordings"
    output_folder_path = repo_root_path / "nwbfiles"
    assets_folder_path = Path(__file__).parent / "assets"

    # Metadata files
    session_metadata_path = assets_folder_path / "session_metadata_cn_amended.csv"

    # List of all subjects
    subject_list = ["#1", "#2", "#3", "#4", "#5", "#6"]

    for selected_subject in tqdm(subject_list, desc="Converting subjects"):
        convert_subject(
            subject_id=selected_subject,
            session_metadata_path=session_metadata_path,
            data_folder_path=data_folder_path,
            output_folder_path=output_folder_path,
        )


if __name__ == "__main__":
    main()
