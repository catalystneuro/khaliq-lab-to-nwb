# Conversion Notes - Khaliq Lab ABF to NWB

## Project Overview
- **Start Date:** March 2025 (data received)
- **Budget:** 250 hours total (50 hours used by Luiz, ~200 hours remaining)
- **Data Status:** Under embargo - NOT for public release
- **Primary Contact:** Lorenzo Sansalone (left NIH June 13, 2025), now Zayd Khaliq and Alejandra Boronat-Garcia

## Data Files

### ABF Files (Primary Focus)
Only ABF files should be processed. Bruker files were not delivered correctly (student responsible left the lab).

#### File Naming Convention (Protocol Prefixes)
- **SF** - Spontaneous Firing
- **CS** - Current Steps
- **ADP** - After Depolarization (afterdepolarization)

#### Critical Issue: File Naming
- **DO NOT RENAME ABF FILES** - Renaming corrupts the files (discovered March 25, 2025)
- Always use original file names from acquisition
- If files are renamed, they become unreadable/corrupted
- This issue occurred multiple times during data sharing (March 25, May 31, June 12)

### Data Corruption
- Do not try to read corrupted ABF files
- Just annotate that there are corrupted parts using `invalid_times` table
- Some files have blocks of corrupted data (high amplitude noise) in specific sweeps
- Corruption identified in initial samples:
  - L_S1C3/ADP_21610017.abf – 2 sweeps
  - L_S1C3/SF_21610014.abf – 1 sweep
  - L_S4C1/ADP_21610043.abf – 1 sweep
  - L_S4C1/SF_21610041.abf – 7 sweeps

### Reading Tools
- Primary: pyABF (Python package for reading ABF files)
- Alternative verification: Clampfit, IgorPro (used by lab)

## Dataset Scope

### Expected Data
**Initial Scope (March 24, 2025):**
- Electrophysiology: Hundreds of cells expected
- Calcium imaging: Tens of cells expected
- Pharmacology data
- Subject metadata

**Actual Data (Updated March 24, 2025):**
- **Electrophysiology:** ~50 cells (no pharmacology)
- **Calcium Imaging:** ~12 cells (with pharmacology)
- Subject metadata and TTL pulses: Available
- **Note:** Not all recordings are sufficient quality for final analysis

### Data Organization
- Experiments are NOT standardized - differences between cells, protocols, and subjects
- Data collected gradually from multiple lab colleagues
- Quality filtering applied before sharing

## Metadata Requirements

### Session Metadata
Must include for each recording:
- Cell unique ID
- Recording date
- Neuron type (e.g., Dopaminergic)
- Anatomical region (e.g., Substantia Nigra)
- Protocol type (SF, CS, ADP)

### Subject Metadata
Required fields:
- Subject unique ID
- Species (e.g., Macaca mulatta)
- Sex
- Age (in years, converted to days for NWB)

### Cell Features
- **AIS** (Axon Initial Segment) - Specific branching of the soma that turns into an axon; this is a cell feature

## Data Sharing

### Platform
- **Globus** file transfer system
- Endpoint: https://app.globus.org/file-manager?origin_id=a0e2af7a-11c6-4118-b15a-f53c219e05e1&origin_path=%2FMy+Drive%2Fdata%2FKhaliq-CN-data-share%2F
- Instructions: https://catalystneuro.com/guides/globus-guide.md

### Access
- Initial: Lorenzo Sansalone
- Current: Zayd Khaliq (khaliqzm@nih.gov) - write access granted July 11, 2025

## Timeline

### Key Dates
- **November 15, 2024:** Initial scope of work sent
- **March 14, 2025:** Data sharing approved by NIH leadership
- **March 17, 2025:** First 6 ABF files uploaded
- **March 25, 2025:** File renaming corruption issue discovered and resolved
- **May 13, 2025:** Most electrophysiology recordings uploaded with Excel metadata
- **May 31, 2025:** Second round of corruption issues identified
- **June 13, 2025:** Lorenzo's last day at NIH
- **July 2025:** Zayd taking over data sharing responsibilities

## Outstanding Items
- Calcium imaging data still being organized (as of May 13, 2025)
- Ongoing data collection for paper in progress
- Alejandra Boronat-Garcia collecting anatomical/transcriptomic data

## Notes
- No publication associated with dataset yet (as of March 21, 2025)
- Data is partial - only what will likely go into future manuscript
- Paper submission in progress (as of July 2, 2025)

## Metadata Additions

### Automatically Extracted from File Structure

The following metadata is automatically extracted from the file path and added to NWB files:

#### Protocol Type (from parent directory)
- **SF**: Spontaneous Firing
- **CS**: Current Steps
- **ADP**: After Depolarization (afterdepolarization)

**Added to NWB:**
- `NWBFile.protocol`: Full name with abbreviation (e.g., "Spontaneous Firing (SF)")
- `NWBFile.keywords`: Protocol abbreviation included

#### AIS Component Status (from directory structure)
- **"with AIS component"**: Cells with Axon Initial Segment component
- **"without AIS component"**: Cells without AIS component

**Added to NWB:**
- `NWBFile.keywords`: AIS status included
- `NWBFile.session_description`: Includes AIS status
- `Subject.description`: Includes AIS status
- `IntracellularElectrode.location`: Includes AIS status

### Metadata Fields Populated

#### NWBFile
- `session_start_time`: From CSV (Date column), timezone: America/New_York (EST/EDT)
- `session_description`: "Intracellular recording from [Neuron Type] neuron ([AIS status]) in [Anatomical Region]"
- `protocol`: "[Protocol Full Name] ([Protocol Abbreviation])"
- `keywords`: [Protocol, AIS status, Neuron Type, Anatomical Region, "patch-clamp", "intracellular"]
- `experimenter`: ["Zayd Khaliq"]
- `lab`: "Cellular Neurophysiology Section"
- `institution`: "National Institute of Neurological Disorders and Stroke, National Institutes of Health"
- `experiment_description`: "Intracellular electrophysiology recording from [Neuron Type] neurons"

#### Subject
- `subject_id`: From CSV (Animal ID column)
- `species`: From CSV, mapped to Latin name (e.g., "Macaca mulatta")
- `sex`: From CSV
- `age`: From CSV, converted to ISO 8601 format (e.g., "P3577D")
- `description`: "[Neuron Type] neuron [AIS status]"

#### IntracellularElectrode
- `location`: "[Anatomical Region] ([AIS status])" (e.g., "Substantia Nigra (with AIS component)")
- `description`: "Intracellular electrode for [Neuron Type] neuron recording"
- `cell_id`: "[Neuron Type]_###" (e.g., "Dopaminergic_001", "Dopaminergic_042")

#### Invalid Times (if file is corrupted)
- Custom columns: "reason", "severity"
- Placeholder intervals for files marked as corrupted
