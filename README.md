# khaliq-lab-to-nwb
NWB conversion scripts for Khaliq lab data to the [Neurodata Without Borders](https://nwb-overview.readthedocs.io/) data format.

**Requirements:** Python 3.13+

## Installation
The package can be installed directly from GitHub, which has the advantage that the source code can be modified if you need to amend some of the code we originally provided to adapt to future experimental differences.
To install the conversion from GitHub you will need to use `git` ([installation instructions](https://github.com/git-guides/install-git)). We also recommend the installation of `conda` ([installation instructions](https://docs.conda.io/en/latest/miniconda.html)) as it contains
all the required machinery in a single and simple install.

From a terminal (note that conda should install one in your system) you can do the following:

```
git clone https://github.com/catalystneuro/khaliq-lab-to-nwb
cd khaliq-lab-to-nwb
conda env create --file make_env.yml
conda activate khaliq_lab
```

This creates a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) which isolates the conversion code from your system libraries.
We recommend that you run all your conversion related tasks and analysis from the created environment in order to minimize issues related to package dependencies.

Alternatively, if you want to avoid conda altogether (for example if you use another virtual environment tool)
you can install the repository with the following commands using only pip:

```
git clone https://github.com/catalystneuro/khaliq-lab-to-nwb
cd khaliq-lab-to-nwb
pip install -e .
```

Or using [uv](https://docs.astral.sh/uv/) (recommended for faster dependency resolution):

```
git clone https://github.com/catalystneuro/khaliq-lab-to-nwb
cd khaliq-lab-to-nwb
uv sync
```

Note:
both pip and uv methods install the repository in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs).

## Repository structure
Each conversion is organized in a directory of its own in the `src` directory:

    khaliq-lab-to-nwb/
    ├── LICENSE
    ├── make_env.yml
    ├── pyproject.toml
    ├── README.md
    ├── notebooks/
    │   └── 001693_demo.ipynb
    └── src/
        └── khaliq_lab_to_nwb/
            └── conversion_abf/
                ├── assets/
                ├── convert_session.py
                ├── icephys_neo_interface.py
                ├── protocol_epochs.py
                ├── saturation_annotation.py
                ├── detect_saturation.py
                └── conversion_notes.md

### conversion_abf

For the `conversion_abf` conversion, you can find a directory located in `src/khaliq_lab_to_nwb/conversion_abf`. This conversion handles intracellular current-clamp recordings (ABF format) from rhesus macaque substantia nigra pars compacta dopaminergic neurons, published as [DANDI:001693](https://dandiarchive.org/dandiset/001693). Inside the conversion directory you can find the following files:

* `convert_session.py`: main conversion script (one NWB file per `(subject, recording_date)`)
* `icephys_neo_interface.py`: Neo `AxonRawIO` adapter that builds `CurrentClampSeries` and `CurrentClampStimulusSeries`
* `protocol_epochs.py`: SF/CS/ADP epoch annotations
* `saturation_annotation.py`: voltage-saturation `invalid_times` flagging
* `detect_saturation.py`: standalone saturation detector and plotting CLI
* `conversion_notes.md`: notes on protocols, scientific context, and the conversion process
* `assets/`: session metadata CSVs

The conversion [notes](src/khaliq_lab_to_nwb/conversion_abf/conversion_notes.md) contain information about the protocols (SF, CS, ADP), the AIS partitioning, and the corruption-tracking workflow.

### Running a specific conversion

You can run the conversion with the following command:

```
python src/khaliq_lab_to_nwb/conversion_abf/convert_session.py
```

## NWB tutorials

A Jupyter notebook demonstrates how to use the NWB files created by the conversion script:

* **conversion_abf**: `notebooks/001693_demo.ipynb`

The notebook streams a session from DANDI by default (DANDI auth required while embargoed) and includes a local-file fallback as commented code in the load cell. It is also under review as a PR against the DANDI example-notebooks repository: [dandi/example-notebooks#148](https://github.com/dandi/example-notebooks/pull/148).

You might need to install `jupyter` before running the notebooks:

```
pip install jupyter
jupyter lab
```
