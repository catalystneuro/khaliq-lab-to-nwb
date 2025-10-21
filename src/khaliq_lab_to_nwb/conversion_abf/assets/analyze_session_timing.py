"""
Analyze session start times for all ABF files to determine temporal grouping.

This script scans all ABF files in the data directory and extracts:
- AIS component status (with/without AIS)
- Recording date folder
- Protocol type (SF, CS, ADP)
- File path
- Session start time

The output is saved as a CSV file for analysis of whether files from the same
date should be grouped into a single experimental session or kept separate.
"""

import csv
from pathlib import Path

from neo.rawio import AxonRawIO


def extract_ais_status(file_path: Path) -> bool:
    """Extract whether the cell has an AIS (Axon Initial Segment) component."""
    path_str = str(file_path)

    if "with AIS component" in path_str:
        return True
    elif "without AIS component" in path_str:
        return False
    else:
        raise ValueError(f"Cannot determine AIS status from path: {file_path}")


def extract_date_folder(file_path: Path) -> str:
    """Extract the date folder name from the file path."""
    # Navigate up from file to find the date folder
    # Structure: .../Cells with AIS component/DATE/CELL/PROTOCOL/file.abf
    parts = file_path.parts

    # Find "Cells with AIS component" or "Cells without AIS component"
    for i, part in enumerate(parts):
        if "Cells with" in part or "Cells without" in part:
            if i + 1 < len(parts):
                return parts[i + 1]

    return "unknown"


def extract_cell_folder(file_path: Path) -> str:
    """Extract the cell folder name (e.g., C1, C12) from the file path."""
    # Navigate up from file to find the cell folder
    # Structure: .../DATE/CELL/PROTOCOL/file.abf
    parts = file_path.parts

    # The cell folder is typically 2 levels up from the file
    if len(parts) >= 3:
        return parts[-3]

    return "unknown"


def analyze_all_sessions(data_folder: Path, output_csv: Path) -> None:
    """
    Analyze all ABF files and save timing information to CSV.

    Parameters
    ----------
    data_folder : Path
        Root folder containing electrophysiology recordings.
    output_csv : Path
        Path where the output CSV file will be saved.
    """
    # Find all ABF files
    abf_files = sorted(data_folder.glob("**/*.abf"))

    print(f"Found {len(abf_files)} ABF files to analyze")

    # Collect data for all files
    rows = []

    for i, abf_file in enumerate(abf_files, 1):
        if i % 10 == 0:
            print(f"Processing file {i}/{len(abf_files)}...")

        try:
            # Extract metadata from path
            ais_status = extract_ais_status(abf_file)
            date_folder = extract_date_folder(abf_file)
            cell_folder = extract_cell_folder(abf_file)
            protocol = abf_file.parent.name

            # Read session start time from ABF file
            reader = AxonRawIO(filename=str(abf_file))
            reader.parse_header()
            session_start_time = reader.raw_annotations["blocks"][0]["rec_datetime"]

            # Format datetime to remove microseconds (only show up to seconds)
            session_start_time_str = session_start_time.strftime("%Y-%m-%dT%H:%M:%S")

            # Store row data
            rows.append(
                {
                    "ais_status": ais_status,
                    "date_folder": date_folder,
                    "cell_folder": cell_folder,
                    "protocol": protocol,
                    "file_name": abf_file.name,
                    "session_start_time": session_start_time_str,
                }
            )

        except Exception as e:
            print(f"Error processing {abf_file}: {e}")
            continue

    # Sort by AIS status, date, cell, and session start time
    rows.sort(
        key=lambda x: (
            x["ais_status"],
            x["date_folder"],
            x["cell_folder"],
            x["session_start_time"],
        )
    )

    # Write to CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(output_csv, "w", newline="") as f:
        fieldnames = [
            "ais_status",
            "date_folder",
            "cell_folder",
            "protocol",
            "file_name",
            "session_start_time",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(rows)

    print(f"\nAnalysis complete! Results saved to: {output_csv}")
    print(f"Total files processed: {len(rows)}")


def main():
    """Main entry point for the script."""
    # Define paths
    repo_root = Path(__file__).parent.parent.parent.parent.parent
    data_folder = repo_root / "data" / "Electrophysiology recordings"
    output_csv = Path(__file__).parent / "session_timing_analysis.csv"

    if not data_folder.exists():
        print(f"Error: Data folder not found at {data_folder}")
        return

    analyze_all_sessions(data_folder, output_csv)


if __name__ == "__main__":
    main()
