# affecttracker_validation

`[Last update: September 15, 2025]`

***
    Period:     2023-04 - 2025-09
    Status:     published

    Author(s):  Antonin Fourcade
    Contact:    antonin.fourcade@maxplanckschools.de

***

## Project description

Scripts for the validation of the AffectTracker, a tool to continously rate 2D (valence, arousal) affective experiences in real-time using Unity

Fourcade, A., Malandrone, F., Roellecke, L., Ciston, A., de Mooij, J., Villringer, A., Carletto S., Gaebler, M. (2024, December 16). AffectTracker: Real-time continuous rating of affective experience in immersive Virtual Reality. PsyArXiv preprint https://doi.org/10.31234/osf.io/xemwb

https://github.com/afourcade/AffectTracker

Data repository: https://doi.org/10.17617/3.QPNSJA

## Project structure

```
affecttracker_validation/
├── README.md                           # Project overview and setup instructions
├── requirements.txt                    # Python dependencies
├── .gitignore                         # Version control exclusions
│
├── code/                              # All source code
│   ├── README.md                      # Code-specific documentation
│   ├── configs/                       # Configuration files
│   │   ├── config.toml               # Main configuration parameters
│   │   └── private_config.toml       # Local/sensitive overrides
│   │
│   └── affecttracker_validation/      # Main Python package
│       ├── __init__.py
│       ├── configs.py                 # Configuration loader/manager
│       │
│       ├── comparison_select_eval/    # Comparison Selection vs. Evaluation studies
│       │   ├── __init__.py
│       │   ├── comparison_phase_survey.py    # Survey phase comparisons
│       │   └── comparison_phase_variability.py # Variability phase comparisons
│       │
│       ├── evaluation/                # Code for Evaluation study
│       │   ├── __init__.py
│       │   │
│       │   ├── datavisualization/    # Plotting and visualization
│       │   │   ├── __init__.py
│       │   │   ├── cr_plot.py        # Continuous ratings plots
│       │   │   ├── cri_sr_corr_plot.py # CRi-SR correlation plots
│       │   │   ├── motion_plot.py    # Motion data visualization
│       │   │   ├── sr_plot.py        # Summary ratings plots
│       │   │   ├── survey_plot.py    # Survey response plots
│       │   │   └── survey_plot_extras.py # Additional survey visualizations
│       │   │
│       │   ├── modeling/              # Statistical analysis & modeling
│       │   │   ├── __init__.py
│       │   │   ├── control_analyses.py     # Control analyses
│       │   │   ├── corr_cr_sr.py           # CRi-SR correlations
│       │   │   ├── cr_freq_change.py       # CR change frequency analysis
│       │   │   └── descriptive_stats.py    # Summary statistics
│       │   │
│       │   └── preprocessing/         # Data cleaning & preparation
│       │       ├── __init__.py
│       │       ├── bids_formatting.py       # BIDS format conversion
│       │       ├── check_length_sequence_videos.py # Video validation
│       │       ├── cri_sr.py               # CRi & SR data processing
│       │       ├── import_cr.py            # Continuous ratings import
│       │       ├── import_motion.py        # Motion data import
│       │       ├── import_PRE_survey.py    # Pre-experiment survey import
│       │       ├── import_POST_survey.py   # Post-experiment survey import
│       │       ├── preprocess_PRE_survey.py # Pre-survey preprocessing
│       │       └── preprocess_POST_survey.py # Post-survey preprocessing
│       │
│       └── selection/                 # Code for Selection study
│           ├── __init__.py
│           │
│           ├── datavisualization/     # Plotting and visualization
│           │   ├── __init__.py
│           │   ├── cr_plot.R          # R-based CR plotting
│           │   ├── cr_plot_extra.R    # extra plots
│           │   └── radar_plot.py      # Radar/spider plots
│           │
│           ├── modeling/              # Statistical analysis & modeling
│           │   ├── __init__.py
│           │   ├── assessment_others.R     # Survey analysis (R)
│           │   ├── corr_cr_sr.R            # Correlations CR-SR (R)
│           │   ├── descriptive_stats.py    # Descriptive stats
│           │   ├── invasiveness.R          # Invasiveness analysis (R)
│           │   └── util.R                  # R utilities
│           │
│           └── preprocessing/         # Data cleaning & preparation
│               ├── __init__.py
│               ├── cri_rs.py               # CRi & SR data processing
│               ├── import_cr.py            # Continuous ratings import
│               ├── import_questionnaire.py # Survey import
│               ├── preprocess_questionnaire.py # Survey preprocessing


```
## Install research code as package

In case, there is no project-related virtual / conda environment yet, create one for the project:

```shell
conda create -n affe_3.11 python=3.11
```

And activate it:

```shell
conda activate affe_3.11
```

Then install the code of the research project as python package:

```shell
# assuming your current working directory is the project root
pip install -e ".[develop]"
```

**Note**: The `-e` flag installs the package in editable mode,
i.e., changes to the code will be directly reflected in the installed package.
Moreover, the code keeps its access to the research data in the underlying folder structure.
Thus, the `-e` flag is recommended to use.

## Publications
Please cite:

Fourcade, A., Malandrone, F., Roellecke, L., Ciston, A., de Mooij, J., Villringer, A., Carletto S., Gaebler, M. (2025). AffectTracker: real-time continuous rating of affective experience in immersive virtual reality. doi: 10.3389/frvir.2025.1567854

## Contributors/Collaborators

Antonin Fourcade
Lucy Roellecke
