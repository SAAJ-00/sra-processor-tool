import os
from pathlib import Path

DEFAULT_CONFIG = {
    'threads': 2,
    'max_size': '30G',
    'output_dir': Path('.').absolute(),
    'keep_temp': False,
    'trim_params': {
        'quality_phred': 20,
        'min_length': 50,
        'cut_window_size': 4,
        'cut_mean_quality': 20,
        'disable_adapter_trimming': False,
        'disable_quality_filtering': False,
        'disable_length_filtering': False
    }
}