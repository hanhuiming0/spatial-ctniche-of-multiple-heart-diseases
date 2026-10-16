import os
from pathlib import Path

ST2_ROOT = Path(os.environ.get('ST2_ROOT', '/path/to/ST2'))
HEART_STUDY_ROOT = Path(os.environ.get('HEART_STUDY_ROOT', '/path/to/Heart_Study'))
DATASHARE_ROOT = Path(os.environ.get('DATASHARE_ROOT', '/path/to/shared_data'))
HMM_DATA_ROOT = Path(os.environ.get('HMM_DATA_ROOT', '/path/to/reference_data'))
PROJECT_ROOT = Path(os.environ.get('PROJECT_ROOT', '.'))
LEGACY_RESULTS_ROOT = Path(os.environ.get('LEGACY_RESULTS_ROOT', '.'))
EXTERNAL_DATA_ROOT = Path(os.environ.get('EXTERNAL_DATA_ROOT', '.'))
REFERENCE_GENESET_ROOT = Path(os.environ.get('REFERENCE_GENESET_ROOT', '.'))
FIGURE_OUTPUT_ROOT = Path(os.environ.get('FIGURE_OUTPUT_ROOT', PROJECT_ROOT / 'outputs'))
