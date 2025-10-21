# Conversion Notes - Khaliq Lab ABF to NWB

## Project Overview
- **Start Date:** March 2025 (data received)
- **Budget:** 250 hours total (50 hours used by Luiz, ~200 hours remaining)
- **Data Status:** Under embargo - NOT for public release
- **Primary Contact:** Lorenzo Sansalone (left NIH June 13, 2025), now Zayd Khaliq and Alejandra Boronat-Garcia

## Scientific Context

### Research Focus
This dataset contains intracellular electrophysiology recordings from dopaminergic neurons in the substantia nigra, examining intrinsic firing properties and responses to current injection protocols (Spontaneous Firing, Current Steps, After Depolarization).

### Axon Initial Segment (AIS)

The **Axon Initial Segment** is the specialized region where the axon emerges from the soma and action potentials are initiated. In substantia nigra dopaminergic neurons, AIS properties show significant cell-to-cell variability.

**Key Characteristics:**
- **Structure:** Variable length (12-60 μm) and diameter (0.2-0.8 μm) with ~50% tapering
- **Function:** The tapering diameter defines a trigger zone for action potentials
- **Molecular:** Predominantly expresses Nav1.2 sodium channels (vs Nav1.6 in other neurons)
- **Firing:** Nav1.2 expression explains the characteristically high threshold for firing and low discharge rate

**Functional Implications:**
Despite large AIS morphological variations (length, position), firing properties remain remarkably consistent across neurons. This robustness is explained by strong somatodendritic excitability that dominates over AIS contributions to neuronal output.

### Relevant Publications

**Publications by Lab Members (Khaliq, Sansalone) - Methodological References:**

1. **Cobb-Lewis DE, Sansalone L, Khaliq ZM** (2023). Contributions of the Sodium Leak Channel NALCN to Pacemaking of Medial Ventral Tegmental Area and Substantia Nigra Dopaminergic Neurons. *J Neurosci* 43(41):6841-6853. DOI: [10.1523/JNEUROSCI.0930-22.2023](https://doi.org/10.1523/JNEUROSCI.0930-22.2023) [PMC10573758]
   - **Relevance: HIGH** - Most similar protocols to current dataset
   - **Protocols Used:** Spontaneous firing (cell-attached), current injection (15s hyperpolarizing ramps)
   - **Equipment:** MultiClamp 700B amplifier, Digidata 1590 digitizer (Molecular Devices)
   - **Electrodes:** Single electrode, 3-7 MΩ borosilicate pipettes
   - **Temperature:** 34°C
   - **Data Acquisition:** 20 kHz sampling, 10 kHz filtering

2. **Evans RC, Twedell EL, Zhu M, Ascencio J, Zhang R, Khaliq ZM** (2020). Functional Dissection of Basal Ganglia Inhibitory Inputs onto Substantia Nigra Dopaminergic Neurons. *Cell Rep* 32(11):108156. DOI: [10.1016/j.celrep.2020.108156](https://doi.org/10.1016/j.celrep.2020.108156) [PMC9887718]
   - **Relevance: MEDIUM-HIGH** - Contains ADP classification and similar protocols
   - **Cell Classification:** "ADP (rebound-ready)" vs "non-ADP (non-rebounding)" based on afterdepolarization
   - **Protocols Used:** Spontaneous firing, current-clamp, rebound firing after hyperpolarization
   - **Electrodes:** Single electrode, 2-7 MΩ borosilicate pipettes
   - **Temperature:** 31-34°C
   - **Note:** ADP here refers to rebound afterdepolarization, not axon initial segment

3. **Evans RC, Zhu M, Khaliq ZM** (2017). Dopamine Inhibition Differentially Controls Excitability of Substantia Nigra Dopamine Neuron Subpopulations through T-Type Calcium Channels. *J Neurosci* 37(13):3704-3720. DOI: [10.1523/JNEUROSCI.0117-17.2017](https://doi.org/10.1523/JNEUROSCI.0117-17.2017) [PMC5373143]
   - **Relevance: HIGH** - Detailed afterdepolarization protocols
   - **Protocols Used:** Spontaneous firing, afterdepolarization (20 Hz AP triplet after hyperpolarization), rebound firing
   - **Electrodes:** Single electrode, 2-7 MΩ borosilicate pipettes
   - **Equipment:** Two-photon imaging (Bruker microscope, Mai Tai Ti:sapphire laser at 810 nm)
   - **Temperature:** Room temp (~22°C) or 30-32°C
   - **Analysis Software:** Igor (Wavemetrics)

4. **Sansalone L, Twedell EL, Evans RC, Boronat-Garcia A, Zhang R, Khaliq ZM** (2025). Corticonigral projections recruit substantia nigra pars lateralis dopaminergic neurons for auditory threat memories. *Nat Commununications* 16(1):8384. DOI: [10.1038/s41467-025-63132-8](https://doi.org/10.1038/s41467-025-63132-8) [PMC11580856]
   - **Relevance: LOW** - Different focus (projection mapping vs intrinsic properties)
   - **Protocols Used:** Spontaneous firing, current injections, afterhyperpolarization
   - **Note:** Focus on SNL vs SNc differences, not AIS categorization

5. **Ding S, Wei W, Zhou FM** (2011). Subcellular Patch-clamp Recordings from the Somatodendritic Domain of Nigral Dopamine Neurons. *J Vis Exp* (49):2769. DOI: [10.3791/2769](https://doi.org/10.3791/2769) [PMC5226116]
   - **Relevance: MEDIUM** - Detailed dual electrode protocols
   - **Equipment:** Axopatch 200B amplifiers (dual recordings possible)
   - **Electrodes:** Soma (6-10 MΩ), Dendrite (8-19 MΩ)
   - **Protocols:** Current steps (depolarizing/hyperpolarizing), dual whole-cell recordings
   - **Note:** Specialized dendritic recording techniques

**AIS Background Literature:**

6. **González-Cabrera C et al.** (2017). Characterization of the axon initial segment of mice substantia nigra dopaminergic neurons. *J Comp Neurol* 525(16):3529-3542. DOI: [10.1002/cne.24288](https://doi.org/10.1002/cne.24288)
   - **Key Finding:** AIS shows variable length (12-60 μm) and diameter (0.2-0.8 μm)
   - **Molecular:** Predominantly Nav1.2 sodium channels (explains high firing threshold)

7. **Moubarak E et al.** (2019). Robustness to axon initial segment variation is explained by somatodendritic excitability in rat substantia nigra dopaminergic neurons. *J Neurosci* 39(26):5044-5065. DOI: [10.1523/JNEUROSCI.2697-18.2019](https://doi.org/10.1523/JNEUROSCI.2697-18.2019)
   - **Key Finding:** Despite large AIS morphological variations, firing properties remain consistent
   - **Conclusion:** Somatodendritic excitability dominates over AIS contributions

**Dataset Context:**
The recordings compare cells "with AIS component" vs "without AIS component". This binary classification does NOT appear in any published literature from the Khaliq lab or AIS research literature, which describes AIS as having **variable morphology** (size, length, position) rather than being present/absent. This suggests the current dataset represents **unpublished work** investigating a novel aspect of AIS structure or function in dopaminergic neurons.

### Typical Recording Protocols from Khaliq Lab Publications

Based on published methods from papers 1-5 above, typical electrophysiology protocols include:

**Protocol Descriptions:**

1. **Spontaneous Firing (SF):**
   - **Mode:** Cell-attached or whole-cell current-clamp
   - **Duration:** 30 seconds to several minutes
   - **Purpose:** Measure baseline pacemaking activity, firing rate, regularity (CV ISI)
   - **Temperature:** 30-34°C (physiological) or room temperature
   - **Analysis:** Firing frequency, interspike interval statistics

2. **Current Steps (CS):**
   - **Mode:** Whole-cell current-clamp
   - **Protocol:** Hyperpolarizing and depolarizing current injections
   - **Duration:** Typically 1-15 seconds per step
   - **Amplitude:** Range from -165 pA to +200 pA in increments (e.g., 20-80 pA steps)
   - **Purpose:** Generate frequency-current (f-I) curves, measure gain (slope) and rheobase (threshold)
   - **Variants:**
     - Short steps (100-500 ms) for action potential characterization
     - Long ramps (15 s) for subthreshold excitability
   - **Analysis:** Input resistance, firing frequency vs current, maximal firing rate

3. **After Depolarization (ADP):**
   - **Mode:** Whole-cell current-clamp
   - **Protocol:** Hyperpolarization followed by current release or depolarizing pulse
   - **Hyperpolarization duration:** 0.9-1 second
   - **Voltage range:** Varied from -90 to -50 mV
   - **Purpose:** Test for rebound excitation, low-threshold calcium currents (T-type), AHP characterization
   - **Trigger protocol:** "20 Hz AP triplet" after hyperpolarization (Evans 2017)
   - **Analysis:**
     - Presence/absence of afterdepolarization
     - Rebound delay (time to first spike)
     - Dendritic Ca²⁺ signals
     - Classification as "ADP (rebound-ready)" vs "non-ADP"

**Standard Equipment Configuration:**

- **Amplifiers:** MultiClamp 700B (Molecular Devices) or Axopatch 200B (dual recordings)
- **Digitizers:** Digidata 1590 (Molecular Devices)
- **Sampling Rate:** 20 kHz
- **Filtering:** 10 kHz low-pass
- **Microscope:** Olympus BX51WI or BX61W1 upright with IR-DIC/DGC optics
- **Camera:** Scientific CMOS (pco.edge 4.2) or similar
- **Objective:** 40× water immersion
- **Temperature Control:** In-line heater maintaining 30-34°C
- **Perfusion:** ~2 ml/min oxygenated ACSF

**Standard Solutions:**

- **Internal (patch pipette):**
  - Base: KMeSO₃ (120-122 mM), NaCl (9 mM), MgCl₂ (1.8-2 mM)
  - ATP/GTP: Mg-ATP (2-4 mM), Na-GTP (0.3-0.5 mM)
  - Buffering: HEPES (9-10 mM), EGTA (0.1-0.45 mM)
  - Energy: Phosphocreatine (5-14 mM)
  - pH 7.35 (KOH)

- **External (ACSF):**
  - Standard: NaCl (125 mM), NaHCO₃ (25 mM), KCl (2.5-3.5 mM), glucose (10-25 mM)
  - Divalent cations: MgCl₂ (1-2 mM), CaCl₂ (2 mM)
  - Oxygenation: 95% O₂ / 5% CO₂
  - pH 7.4, 300-310 mOsm

**Recording Configuration:**
- **Single electrode** for most protocols (cell-attached or whole-cell)
- **Dual electrode** possible for specialized dendritic recordings (Axopatch 200B × 2)
- **Pipette resistance:** 2-7 MΩ (soma), 8-19 MΩ (dendrite)
- **Access resistance:** <25 MΩ, compensated 30-70%

**Data Acquisition Software:**
- Primary: pCLAMP (Molecular Devices)
- Analysis: Igor Pro (Wavemetrics), Prism (GraphPad), custom routines in NeuroMatic toolkit

**ABF File Format:**
All Khaliq lab data uses Axon Binary Format (ABF) version 2.x, which stores:
- Recording date/time (uFileStartDate + uFileStartTimeMS)
- Multiple sweeps per protocol
- Channel configurations (voltage, current)
- Protocol metadata (sampling rate, filter settings, gains)

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

## Typical Experimental Session Description

### Overview of a Complete Recording Session

Based on the dataset structure and published Khaliq lab methods, a typical experimental session follows this workflow:

**Session Timeline (typical duration: 30-60 minutes per cell):**

1. **Preparation Phase (Before recording)**
   - Acute brain slice preparation from Rhesus macaque (Macaca mulatta)
   - Target: Substantia nigra pars compacta (SNc) dopaminergic neurons
   - Slice incubation at 34°C in oxygenated ACSF for 30-60 minutes recovery
   - Tissue classification: Cell identified as "with AIS component" or "without AIS component" (classification method unknown - likely based on prior anatomical analysis or functional criteria)

2. **Cell Selection and Targeting (5-10 minutes)**
   - Recording chamber mounted on Olympus BX51WI upright microscope
   - Continuous perfusion with oxygenated ACSF at ~2 ml/min, 30-34°C
   - Visual identification using IR-DIC or Dodt gradient contrast optics (40× objective)
   - Cell selection criteria:
     - Healthy appearance (smooth membrane, visible nucleus)
     - Soma located 10-30 μm below slice surface
     - Long, visible dendrites indicating viable neuron
   - Initial electrophysiological verification:
     - Firing frequency <5 Hz (characteristic of DA neurons)
     - Presence of HCN-mediated sag (Ih current)

3. **Seal Formation and Whole-Cell Access (5-10 minutes)**
   - Borosilicate pipette (3-7 MΩ resistance) filled with K-based internal solution
   - Approach under positive pressure, then gentle seal formation (>1 GΩ)
   - Whole-cell configuration achieved by brief suction
   - Quality control:
     - Series resistance <25 MΩ (compensated 30-70%)
     - Stable resting membrane potential
     - Healthy action potential amplitude and shape

4. **Protocol Sequence (20-40 minutes total)**

   **Protocol 1: Spontaneous Firing (SF) - 2-5 minutes**
   - **File naming:** Prefix "SF" or numeric code (e.g., 21610014.abf)
   - **Recording mode:** Whole-cell current-clamp (I=0 mode)
   - **Duration:** 30 seconds to several minutes of continuous recording
   - **Temperature:** 30-34°C maintained throughout
   - **Measurements:**
     - Baseline firing rate (typically 1-4 Hz for SNc DA neurons)
     - Regularity of pacemaking (coefficient of variation of ISI)
     - Action potential waveform characteristics
     - Spontaneous subthreshold oscillations
   - **Data saved:** Single ABF file with continuous sweeps

   **Protocol 2: Current Steps (CS) - 10-15 minutes**
   - **File naming:** Prefix "CS" or numeric code (e.g., 21610015.abf)
   - **Recording mode:** Whole-cell current-clamp
   - **Protocol structure:**
     - Series of depolarizing and hyperpolarizing current steps
     - Typical range: -165 pA to +200 pA
     - Step increment: 20-80 pA
     - Step duration: 1-15 seconds per step
     - Alternative: Long hyperpolarizing ramps (15 s duration)
   - **Measurements:**
     - Input resistance (from voltage deflection)
     - Frequency-current (f-I) relationship
     - Gain (slope of f-I curve)
     - Rheobase (minimum current to elicit firing)
     - Maximal firing frequency
     - Sag voltage (Ih activation during hyperpolarization)
   - **Data saved:** Single ABF file with multiple sweeps (one per current level)

   **Protocol 3: After Depolarization (ADP) - 10-20 minutes**
   - **File naming:** Prefix "ADP" or numeric code (e.g., 21610017.abf)
   - **Recording mode:** Whole-cell current-clamp
   - **Protocol structure:**
     - Hold neuron at hyperpolarized potential (-90 to -50 mV) for 0.9-1 second
     - Release current (rebound protocol) OR deliver brief depolarizing pulse
     - Alternative: Evoke "20 Hz AP triplet" after hyperpolarization period
     - Repeat at multiple holding voltages
   - **Measurements:**
     - Presence/absence of rebound depolarization
     - Afterhyperpolarization (AHP) amplitude and duration
     - Low-threshold calcium spike (T-type channel activation)
     - Rebound delay (latency to first spike)
     - Rebound burst characteristics
     - Cell classification: "ADP (rebound-ready)" vs "non-ADP"
   - **Data saved:** Single ABF file with multiple sweeps (one per holding voltage)

5. **Session Completion and Documentation (5 minutes)**
   - Visual confirmation that cell remained healthy throughout
   - In some experiments: Biocytin filling for post-hoc morphological reconstruction
   - Recording of final series resistance and membrane properties
   - Files saved with original acquisition names (NEVER renamed due to corruption risk)
   - Session metadata logged: Date, cell number, animal ID, protocols completed

**Temporal Clustering:**
Analysis shows that all three protocols (SF, CS, ADP) are typically completed within a **5-10 minute window** for the same cell, with 98.9% of recordings within the same cell folder occurring <30 minutes apart. This indicates the three ABF files represent **sequential protocols on the same neuron** during a single experimental session.

**Directory Structure After Session:**
```
Electrophysiology recordings/
└── [Cells with AIS component OR Cells without AIS component]/
    └── [Date folder, e.g., "10 June 2021"]/
        └── [Cell folder, e.g., "C1"]/
            ├── SF/
            │   └── 21610014.abf  (Spontaneous Firing)
            ├── CS/
            │   └── 21610015.abf  (Current Steps)
            └── ADP/
                └── 21610017.abf  (After Depolarization)
```

**Equipment State During Session:**
- **Amplifier:** MultiClamp 700B in current-clamp and voltage-clamp modes
- **Digitizer:** Digidata 1590 sampling at 20 kHz, filtering at 10 kHz
- **Acquisition software:** pCLAMP (Molecular Devices) - generates ABF 2.x format files
- **Temperature:** Maintained at 30-34°C via in-line heater with feedback control
- **Perfusion:** Continuous ACSF flow (95% O₂ / 5% CO₂ bubbled)

**Data Quality Indicators:**
- Stable series resistance throughout session (<25 MΩ, <20% change)
- Consistent action potential amplitude
- Stable resting membrane potential
- No signs of seal deterioration or cell damage

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

## Session Timing Analysis

### Extracting Session Start Times

Session start times are extracted from ABF files using the neo library:

```python
from neo.rawio import AxonRawIO

reader = AxonRawIO(filename='file.abf')
reader.parse_header()
session_start_time = reader.raw_annotations['blocks'][0]['rec_datetime']
```

This is the same approach used by neuroconv's AbfInterface, which accesses the underlying ABF metadata fields:
- **ABF 2.x**: Uses `uFileStartDate` (YYYYMMDD) + `uFileStartTimeMS` (milliseconds since midnight)
- **ABF 1.x**: Falls back to `rec_datetime` (only time-of-day, date defaults to 1900-01-01)

All Khaliq Lab data uses ABF 2.x format.

### Recording Session Grouping

Analysis of 136 ABF files across 41 cells reveals that recordings within each cell are temporally clustered:

**Key Findings:**
- **98.9% of recordings** within each cell are less than 30 minutes apart
- Typical pattern: SF → CS → ADP protocols run sequentially within 5-10 minutes total
- Average time between consecutive protocols: 2-6 minutes

**Recommendation:**
Recordings from the same cell/date folder represent a single experimental session where different protocols (SF, CS, ADP) are applied sequentially to the same cell. Consider:
1. Adding a `session_id` metadata field to group related recordings
2. Combining protocol recordings into a single NWB file per cell
3. At minimum, documenting this relationship in NWB metadata

**Analysis Files:**
- `assets/session_timing_analysis.csv` - Full timing data for all files (columns: ais_status, date_folder, cell_folder, protocol, file_name, session_start_time)
- `assets/analyze_session_timing.py` - Script to regenerate analysis
