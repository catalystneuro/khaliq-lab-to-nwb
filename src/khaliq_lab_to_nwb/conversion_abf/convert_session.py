from pathlib import Path
from typing import Literal
from zoneinfo import ZoneInfo

import pandas as pd
from neuroconv.datainterfaces import AbfInterface
from neuroconv.tools import configure_and_write_nwbfile


def load_metadata_from_csv(
    cell_number: int,
    session_metadata_path: Path,
) -> dict:
    """
    Load raw metadata for a specific cell from the session metadata CSV file.

    This function finds and returns the matching row from the CSV as raw, unprocessed data.

    Parameters
    ----------
    cell_number : int
        The cell number to look up (corresponds to "Cell #" column in CSV).
    session_metadata_path : Path
        Path to the session metadata CSV file.

    Returns
    -------
    dict
        Dictionary containing raw CSV data:
        - "date": Recording date (str or datetime)
        - "cell_number": Cell identifier (int)
        - "neuron_type": Type of neuron (str)
        - "anatomical_region": Brain region (str)
        - "animal_id": Subject identifier (str)
        - "animal_species": Species name (str)
        - "sex": Subject sex (str)
        - "age_years": Subject age in years (float)

    Raises
    ------
    ValueError
        If the specified cell number is not found in the session metadata.
    """
    df = pd.read_csv(session_metadata_path)

    # Find the row for this cell
    cell_row = df[df["Cell #"] == cell_number]

    if cell_row.empty:
        raise ValueError(f"Cell #{cell_number} not found in Master List")

    cell_data = cell_row.iloc[0]

    # Return raw CSV data
    return {
        "date": cell_data["Date"],
        "cell_number": int(cell_data["Cell #"]),
        "neuron_type": cell_data["Neuron Type"],
        "anatomical_region": cell_data["Anatomical Region"],
        "animal_id": cell_data["Animal ID"],
        "animal_species": cell_data["Animal Species"],
        "sex": cell_data["Sex"],
        "age_years": float(cell_data["Animal Age (Years)"]),
    }


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


# ============================================================================
# MAIN CONVERSION WORKFLOW
# ============================================================================


def convert_session(
    abf_file_path: Path | str,
    output_folder_path: Path | str,
    session_metadata: dict,
    cell_number: int,
    ais_status: Literal["with AIS component", "without AIS component"],
) -> Path:
    """
    Convert a single ABF electrophysiology session to NWB format.

    This is the main conversion function that orchestrates the entire workflow:
    1. Validates the ABF file exists
    2. Checks corruption status
    3. Initializes NeuroConv ABF interface
    4. Creates NWB file with enriched metadata
    5. Adds corruption annotations if applicable
    6. Writes the NWB file

    Parameters
    ----------
    abf_file_path : Path or str
        Path to the ABF file to convert
    output_folder_path : Path or str
        Path to folder where NWB file will be saved. Directory will be created
        if it doesn't exist. The output filename will be formatted as
        "cell_{number:03d}_{abf_stem}.nwb".
    session_metadata : dict
        Dictionary containing raw CSV metadata
    cell_number : int
        The cell number for this recording session
    ais_status : Literal["with AIS component", "without AIS component"]
        Whether the cell has an Axon Initial Segment component

    Returns
    -------
    Path
        Path to the created NWB file

    Raises
    ------
    FileNotFoundError
        If the ABF file does not exist

    Notes
    -----
    - The function applies a workaround for a NeuroConv bug by setting
      icephys_experiment_type directly in metadata
    - Corruption annotations use the invalid_times table with custom columns
      for "reason" and "severity"
    """
    abf_file_path = Path(abf_file_path)
    output_folder_path = Path(output_folder_path)

    # Validate input file exists
    if not abf_file_path.exists():
        raise FileNotFoundError(f"ABF file not found: {abf_file_path}")

    # Process raw CSV metadata
    # Extract and format the date
    date_value = session_metadata["date"]
    if isinstance(date_value, pd.Timestamp):
        session_start_time = date_value.to_pydatetime()
    else:
        session_start_time = pd.to_datetime(date_value).to_pydatetime()

    # Add timezone info (NIH Bethesda, MD is in US Eastern Time)
    session_start_time = session_start_time.replace(tzinfo=ZoneInfo("America/New_York"))

    # Convert common name to Latin name
    species_mapping = {
        "Rhesus Macaque": "Macaca mulatta",
        "Rhesus macaque": "Macaca mulatta",
    }
    species_latin = species_mapping.get(
        session_metadata["animal_species"], session_metadata["animal_species"]
    )

    # Convert age from years to days in ISO 8601 format
    age_iso = f"P{session_metadata['age_years'] * 365:.0f}D"

    # Extract protocol type from file path
    protocol_type = abf_file_path.parent.name
    valid_protocols = ["SF", "CS", "ADP"]
    if protocol_type not in valid_protocols:
        raise ValueError(
            f"Unknown protocol type '{protocol_type}' in path: {abf_file_path}"
        )

    # Check corruption status
    is_corrupted = check_corruption_status(abf_file_path)

    # Step 1: Initialize the ABF interface
    interface = AbfInterface(file_paths=[abf_file_path])

    # Step 2: Get and customize metadata
    metadata = interface.get_metadata()

    # Map protocol abbreviations to full names
    protocol_names = {
        "SF": "Spontaneous Firing",
        "CS": "Current Steps",
        "ADP": "After Depolarization",
    }
    protocol_full_name = protocol_names.get(protocol_type, protocol_type)

    # Update NWBFile metadata
    metadata["NWBFile"]["session_start_time"] = session_start_time
    metadata["NWBFile"]["session_description"] = (
        f"Intracellular recording from {session_metadata['neuron_type']} neuron "
        f"({ais_status}) in {session_metadata['anatomical_region']}"
    )
    metadata["NWBFile"]["protocol"] = f"{protocol_full_name} ({protocol_type})"
    metadata["NWBFile"]["experimenter"] = ["Khaliq, Zayd", "Sansalone, Lorenze"]
    metadata["NWBFile"]["lab"] = "Cellular Neurophysiology Section"
    metadata["NWBFile"]["institution"] = (
        "National Institute of Neurological Disorders and Stroke, National Institutes of Health"
    )
    metadata["NWBFile"]["experiment_description"] = (
        f"Intracellular electrophysiology recording from {session_metadata['neuron_type']} neurons"
    )
    metadata["NWBFile"]["keywords"] = [
        protocol_type,
        ais_status,
        session_metadata["neuron_type"],
        session_metadata["anatomical_region"],
        "patch-clamp",
        "intracellular",
    ]

    # Update Subject metadata
    metadata["Subject"]["subject_id"] = session_metadata["animal_id"]
    metadata["Subject"]["species"] = species_latin
    metadata["Subject"]["sex"] = session_metadata["sex"]
    metadata["Subject"]["age"] = age_iso
    metadata["Subject"]["description"] = (
        f"{session_metadata['neuron_type']} neuron {ais_status}"
    )

    # Update Icephys Electrode metadata with anatomical location and cell ID
    if "Icephys" in metadata and "Electrodes" in metadata["Icephys"]:
        electrode_location = f"{session_metadata['anatomical_region']} ({ais_status})"
        for electrode in metadata["Icephys"]["Electrodes"]:
            electrode["location"] = electrode_location
            electrode["description"] = (
                f"Intracellular electrode for {session_metadata['neuron_type']} neuron recording"
            )
            electrode["cell_id"] = (
                f"{session_metadata['neuron_type']}_{session_metadata['cell_number']:03d}"
            )

    # WORKAROUND for neuroconv bug: Set icephys_experiment_type in metadata
    # The add_to_nwbfile() method ignores the icephys_experiment_type parameter
    # and only reads from metadata["Icephys"]["Sessions"][i]["icephys_experiment_type"]
    # See: neuroconv_bug_report.md
    for session in metadata["Icephys"]["Sessions"]:
        session["icephys_experiment_type"] = "current_clamp"

    # Step 3: Create NWB file with the interface
    nwbfile = interface.create_nwbfile(metadata=metadata)

    # Step 4: Add corruption annotation ONLY if file is corrupted
    if is_corrupted:
        nwbfile.add_invalid_times_column(
            name="reason",
            description="Type of corruption or artifact detected in the data",
        )
        nwbfile.add_invalid_times_column(
            name="severity", description="Severity level of the corruption"
        )

        # Add a placeholder interval indicating the file has corruption
        # NOTE: Specific time intervals should be determined by detailed analysis
        nwbfile.add_invalid_time_interval(
            start_time=0.0,
            stop_time=0.0,  # Placeholder - actual intervals need to be determined
            reason="file_marked_as_corrupted",
            severity="unknown",
        )

    # Step 5: Write to NWB file
    output_folder_path.mkdir(parents=True, exist_ok=True)

    # Generate output filename
    output_filename = f"cell_{cell_number:03d}_{abf_file_path.stem}.nwb"
    nwbfile_path = output_folder_path / output_filename

    # Step 6: Write to NWB file
    configure_and_write_nwbfile(
        nwbfile=nwbfile, nwbfile_path=nwbfile_path, backend="hdf5"
    )

    return nwbfile_path


def main():
    """
    Main entry point for ABF to NWB conversion.

    This function sets up default paths for the Khaliq Lab repository structure
    and converts an example ABF file. Modify the abf_file_path variable to
    convert different files.

    The function demonstrates typical usage of the convert_session function
    with the expected directory layout.
    """

    # ============================================================================
    # CONFIGURATION
    # ============================================================================

    # Base paths
    repo_root_path = Path(__file__).parent.parent.parent.parent
    data_folder_path = repo_root_path / "data" / "Electrophysiology recordings"
    output_folder_path = repo_root_path / "nwbfiles"
    assets_folder_path = Path(__file__).parent / "assets"

    # Metadata files
    session_metadata_path = assets_folder_path / "session_metadata.csv"

    # Example ABF file path - UPDATE THIS
    ais_type = "with AIS component"  # or "without AIS component"
    folder_date = "10 June 2021"
    cell = "C1"
    protocol = "ADP"
    abf_file_path = (
        data_folder_path / f"{ais_type}/{folder_date}/{cell}/{protocol}/21610017.abf"
    )

    # Extract cell number and load metadata
    cell_number = extract_cell_number_from_path(abf_file_path)
    session_metadata = load_metadata_from_csv(cell_number, session_metadata_path)

    # Extract AIS status from file path
    path_str = str(abf_file_path)
    if "with AIS component" in path_str:
        ais_status = "with AIS component"
    elif "without AIS component" in path_str:
        ais_status = "without AIS component"
    else:
        raise ValueError(f"Cannot determine AIS status from path: {abf_file_path}")

    # Convert the example file
    convert_session(
        abf_file_path=abf_file_path,
        output_folder_path=output_folder_path,
        session_metadata=session_metadata,
        cell_number=cell_number,
        ais_status=ais_status,
    )


if __name__ == "__main__":
    main()
