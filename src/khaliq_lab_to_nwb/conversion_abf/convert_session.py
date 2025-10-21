from pathlib import Path
from zoneinfo import ZoneInfo
import pandas as pd

from neuroconv.datainterfaces import AbfInterface
from neuroconv.tools import configure_and_write_nwbfile


def load_metadata_from_csv(
    cell_number: int,
    session_metadata_path: Path,
) -> dict:
    """
    Load metadata for a specific cell from the session metadata CSV file.

    This function reads the session metadata CSV file containing experimental metadata
    for all cells, locates the row corresponding to the specified cell number, and
    extracts relevant information to populate NWB metadata fields.

    Parameters
    ----------
    cell_number : int
        The cell number to look up in the session metadata (e.g., 1, 2, 3...).
        This corresponds to the "Cell #" column in the CSV file.
    session_metadata_path : Path
        Path to the session metadata CSV file containing columns:
        - "Cell #": Cell identifier
        - "Date": Recording date
        - "Neuron Type": Type of neuron (e.g., "Dopaminergic")
        - "Anatomical Region": Brain region (e.g., "Substantia Nigra")
        - "Animal ID": Subject identifier
        - "Animal Species": Species name (e.g., "Rhesus Macaque")
        - "Sex": Subject sex ("M" or "F")
        - "Animal Age (Years)": Subject age in years

    Returns
    -------
    dict
        Dictionary containing NWB-formatted metadata with keys:
        - "session_start_time": datetime with timezone
        - "session_description": str
        - "experimenter": list[str]
        - "lab": str
        - "institution": str
        - "experiment_description": str
        - "subject": dict with subject_id, species, sex, age

    Raises
    ------
    ValueError
        If the specified cell number is not found in the session metadata
    FileNotFoundError
        If the session metadata CSV file does not exist

    Notes
    -----
    - All timestamps are converted to US Central Time (America/Chicago)
    - Species names are mapped from common names to Latin nomenclature
    - Age is converted from years to ISO 8601 duration format (e.g., "P3577D")
    """
    df = pd.read_csv(session_metadata_path)

    # Find the row for this cell
    cell_row = df[df["Cell #"] == cell_number]

    if cell_row.empty:
        raise ValueError(f"Cell #{cell_number} not found in Master List")

    cell_data = cell_row.iloc[0]

    # Extract and format the date
    date_value = cell_data["Date"]
    if isinstance(date_value, pd.Timestamp):
        session_start_time = date_value.to_pydatetime()
    else:
        session_start_time = pd.to_datetime(date_value).to_pydatetime()

    # Add timezone info (assuming US Central Time for the lab)
    session_start_time = session_start_time.replace(tzinfo=ZoneInfo("America/Chicago"))

    # Convert common name to Latin name
    species_mapping = {
        "Rhesus Macaque": "Macaca mulatta",
        "Rhesus macaque": "Macaca mulatta",
    }
    species_latin = species_mapping.get(cell_data["Animal Species"], cell_data["Animal Species"])

    # Build metadata dictionary
    metadata = {
        "session_start_time": session_start_time,
        "session_description": f"Intracellular recording from {cell_data['Neuron Type']} neuron in {cell_data['Anatomical Region']}",
        "experimenter": ["Khaliq Lab"],
        "lab": "Khaliq Lab",
        "institution": "Northwestern University",
        "experiment_description": f"Electrophysiology recording from {cell_data['Neuron Type']} neurons",
        "subject": {
            "subject_id": cell_data["Animal ID"],
            "species": species_latin,
            "sex": cell_data["Sex"],
            "age": f"P{float(cell_data['Animal Age (Years)']) * 365:.0f}D",  # Convert years to days
        },
    }

    return metadata


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
    corruption_csv_path = Path(__file__).parent / "assets" / "file_corruption_status.csv"
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
        Dictionary containing NWB-formatted metadata with keys:
        - "session_start_time": datetime with timezone
        - "session_description": str
        - "experimenter": list[str]
        - "lab": str
        - "institution": str
        - "experiment_description": str
        - "subject": dict with subject_id, species, sex, age
    cell_number : int
        The cell number for this recording session

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

    # Check corruption status
    is_corrupted = check_corruption_status(abf_file_path)

    # Step 1: Initialize the ABF interface
    interface = AbfInterface(file_paths=[abf_file_path])

    # Step 2: Get and customize metadata
    metadata = interface.get_metadata()

    # Update NWBFile metadata
    metadata["NWBFile"]["session_start_time"] = session_metadata["session_start_time"]
    metadata["NWBFile"]["session_description"] = session_metadata["session_description"]
    metadata["NWBFile"]["experimenter"] = session_metadata["experimenter"]
    metadata["NWBFile"]["lab"] = session_metadata["lab"]
    metadata["NWBFile"]["institution"] = session_metadata["institution"]
    metadata["NWBFile"]["experiment_description"] = session_metadata["experiment_description"]

    # Update Subject metadata
    metadata["Subject"].update(session_metadata["subject"])

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
            name="reason", description="Type of corruption or artifact detected in the data"
        )
        nwbfile.add_invalid_times_column(name="severity", description="Severity level of the corruption")

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
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path, backend="hdf5")


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
    data_folder_path = repo_root_path / "data"
    output_folder_path = repo_root_path / "nwbfiles"
    assets_folder_path = Path(__file__).parent / "assets"

    # Metadata files
    session_metadata_path = assets_folder_path / "session_metadata.csv"

    # Example ABF file path - UPDATE THIS
    abf_file_path = data_folder_path / "Electrophysiology recordings/Cells with AIS component/10 June 2021/C1/ADP/21610017.abf"

    # Extract cell number and load metadata
    cell_number = extract_cell_number_from_path(abf_file_path)
    session_metadata = load_metadata_from_csv(cell_number, session_metadata_path)

    # Convert the example file
    convert_session(
        abf_file_path=abf_file_path,
        output_folder_path=output_folder_path,
        session_metadata=session_metadata,
        cell_number=cell_number,
    )


if __name__ == "__main__":
    main()
