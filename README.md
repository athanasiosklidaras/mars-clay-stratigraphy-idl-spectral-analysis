# Mars spectral analysis pipeline (IDL)

IDL code for the spectral analysis pipeline used in Athanasios Klidaras's PhD research on Mars. The code is also archived as `spectral_analysis_pipeline_IDL_code.zip` in the Mendeley Data record below.

The associated dataset describes analysis of map-projected targeted observations from the Compact Reconnaissance Imaging Spectrometer for Mars (CRISM), including products used to investigate compositional clay stratigraphies and weathering on ancient Mars. The dataset and other supporting files remain available from Mendeley Data; the original CRISM products are available separately from the NASA Planetary Data System.

## What the code does

The pipeline generates clay band-center and band-depth maps, a clay detection mask and clay-phase classifications; calculates spectral parameters; and selects regions of interest and extracts spectra. The source file contains further comments describing its methods and history.

## Contents and requirements

- `Klidaras_2026_global_signatures_weathering.pro`: the IDL source file, unchanged from the Mendeley code archive.
- The source describes use in an IDL/ENVI environment with CRISM input cubes. The full dataset and original CRISM observations are not included in this repository.
- The batch loader contains a machine-specific input path, and output files are written under `Saved data paper2/`. These paths need to be adapted for another environment.

The comments in the source show examples for processing a single observation or the 89-observation database. The code is shared as archived; the examples have not been tested as part of creating this repository.

## Credits

Written by Athanasios (Thanos) Klidaras with input from Briony Horgan and contributions from Chandler Albright. See the source file header for additional context and acknowledgements.

## Archive and citation

Klidaras, Athanasios (2026), “Dataset for Global Spectral Signatures of Weathering on the Surface of Ancient Mars”, Mendeley Data, V1. [https://doi.org/10.17632/s3xj33wfb4.1](https://doi.org/10.17632/s3xj33wfb4.1).

The Mendeley record is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Consult that record for the full dataset, documentation, and attribution details.
